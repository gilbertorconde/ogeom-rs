//! Point clouds, and the text formats that carry them: XYZ and PTS.
//!
//! A cloud is points with what a scanner measured at each: a normal, a
//! colour, and any number of named scalars (intensity, confidence, a
//! deviation an inspection computed). [`PointCloud`] holds them side by
//! side, one entry per point in each.
//!
//! XYZ is a row of numbers per point, separated by spaces, tabs, commas or
//! semicolons; PTS is the same with a line before each scan saying how many
//! rows follow. Neither says what its columns are, so a caller may name
//! them ([`XyzColumn`]); otherwise a header line naming them (`# x y z nx ny
//! nz`, `//X,Y,Z,R,G,B`) is honoured, and otherwise the column count decides,
//! as [`read_xyz`] lays out, and a count that leaves the reading open is
//! refused by name. Neither declares a unit; coordinates come back as
//! written.

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::{Point, Vector};
use std::fmt::Write as _;

/// One named value per point.
#[derive(Debug, Clone, PartialEq)]
pub struct Scalar {
    /// The name the file gives it.
    pub name: String,
    /// One value per point, in point order.
    pub values: Vec<f64>,
}

/// Points, with what was measured at each.
///
/// Every attribute present has exactly one entry per point;
/// [`PointCloud::check`] says whether one does not.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PointCloud {
    /// The points.
    pub points: Vec<Point>,
    /// A normal per point, as the file gives it, where it gives one.
    pub normals: Option<Vec<Vector>>,
    /// An RGBA colour per point in `[0, 1]`, where the file gives one.
    pub colours: Option<Vec<[f64; 4]>>,
    /// Named scalars per point, in the order the file gives them.
    pub scalars: Vec<Scalar>,
}

impl PointCloud {
    /// Bare points, with nothing measured at them.
    #[must_use]
    pub fn new(points: Vec<Point>) -> Self {
        Self {
            points,
            ..Self::default()
        }
    }

    /// The number of points.
    #[must_use]
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether there are no points.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// The values of the scalar named `name`.
    #[must_use]
    pub fn scalar(&self, name: &str) -> Option<&[f64]> {
        self.scalars
            .iter()
            .find(|s| s.name == name)
            .map(|s| s.values.as_slice())
    }

    /// Check that every attribute has one entry per point.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) naming
    /// the attribute whose length differs from the point count.
    pub fn check(&self) -> OgeomResult<()> {
        let n = self.points.len();
        if let Some(normals) = &self.normals
            && normals.len() != n
        {
            ogeom_bail!(Dimension, "{} normals for {n} points", normals.len());
        }
        if let Some(colours) = &self.colours
            && colours.len() != n
        {
            ogeom_bail!(Dimension, "{} colours for {n} points", colours.len());
        }
        for s in &self.scalars {
            if s.values.len() != n {
                ogeom_bail!(
                    Dimension,
                    "{} values of {} for {n} points",
                    s.values.len(),
                    s.name
                );
            }
        }
        Ok(())
    }
}

/// What one column of an XYZ or PTS row holds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum XyzColumn {
    /// The point's `x`.
    X,
    /// The point's `y`.
    Y,
    /// The point's `z`.
    Z,
    /// The normal's `x`.
    NormalX,
    /// The normal's `y`.
    NormalY,
    /// The normal's `z`.
    NormalZ,
    /// Red, from 0 to 255.
    Red,
    /// Green, from 0 to 255.
    Green,
    /// Blue, from 0 to 255.
    Blue,
    /// A named scalar.
    Scalar(String),
    /// A column to pass over.
    Skip,
}

impl XyzColumn {
    /// The column a header word names; anything not a coordinate, a normal
    /// or a colour channel is a scalar of that name.
    fn named(word: &str) -> Self {
        match word.to_ascii_lowercase().as_str() {
            "x" => Self::X,
            "y" => Self::Y,
            "z" => Self::Z,
            "nx" => Self::NormalX,
            "ny" => Self::NormalY,
            "nz" => Self::NormalZ,
            "r" | "red" => Self::Red,
            "g" | "green" => Self::Green,
            "b" | "blue" => Self::Blue,
            _ => Self::Scalar(word.to_owned()),
        }
    }

    fn word(&self) -> &str {
        match self {
            Self::X => "x",
            Self::Y => "y",
            Self::Z => "z",
            Self::NormalX => "nx",
            Self::NormalY => "ny",
            Self::NormalZ => "nz",
            Self::Red => "red",
            Self::Green => "green",
            Self::Blue => "blue",
            Self::Scalar(name) => name,
            Self::Skip => "-",
        }
    }
}

fn fields(line: &str) -> impl Iterator<Item = &str> {
    line.split(|c: char| c.is_whitespace() || c == ',' || c == ';')
        .filter(|f| !f.is_empty())
}

/// A header line naming its columns, `x y z` first.
fn header(line: &str) -> Option<Vec<XyzColumn>> {
    let rest = line
        .strip_prefix("//")
        .or_else(|| line.strip_prefix('#'))?
        .trim_start_matches(|c: char| c == '#' || c == '/' || c.is_whitespace());
    let rest = rest.strip_prefix("columns:").map_or(rest, str::trim_start);
    let columns: Vec<XyzColumn> = fields(rest).map(XyzColumn::named).collect();
    columns
        .starts_with(&[XyzColumn::X, XyzColumn::Y, XyzColumn::Z])
        .then_some(columns)
}

/// The rows of a text cloud, column by column.
struct Table {
    columns: Vec<Vec<f64>>,
    header: Option<Vec<XyzColumn>>,
}

/// Read the numeric rows. `counted` is PTS: a line of one integer opens a
/// scan of that many rows.
fn table(text: &str, counted: bool) -> OgeomResult<Table> {
    let mut columns: Vec<Vec<f64>> = Vec::new();
    let mut header_line = None;
    let mut seen_data = false;
    // PTS: the rows the open scan still owes.
    let mut owed: Option<usize> = None;
    for (number, line) in text.lines().enumerate() {
        let line = line.trim();
        let number = number + 1;
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') || line.starts_with("//") {
            if !seen_data && header_line.is_none() {
                header_line = header(line);
            }
            continue;
        }
        let mut row = fields(line).peekable();
        if counted && let Some(only) = row.peek().copied() {
            let mut probe = fields(line);
            probe.next();
            if probe.next().is_none() {
                let Ok(count) = only.parse::<usize>() else {
                    ogeom_bail!(
                        Construction,
                        "line {number} of the PTS holds one field, {only}, which is not a point count"
                    );
                };
                if let Some(left) = owed.filter(|&n| n > 0) {
                    ogeom_bail!(
                        Construction,
                        "line {number} opens a scan while the one before still owes {left} points"
                    );
                }
                owed = Some(count);
                continue;
            }
            match owed.as_mut() {
                None => ogeom_bail!(
                    Construction,
                    "the PTS has no point count before its first point, on line {number}"
                ),
                Some(0) => ogeom_bail!(
                    Construction,
                    "line {number} is a point past the count its scan gave"
                ),
                Some(n) => *n -= 1,
            }
        }
        let mut width = 0;
        for (k, field) in row.enumerate() {
            let Ok(value) = field.parse::<f64>() else {
                ogeom_bail!(
                    Construction,
                    "line {number} carries {field}, which is not a number"
                );
            };
            if !value.is_finite() {
                ogeom_bail!(
                    Construction,
                    "line {number} carries {field}, which is not finite"
                );
            }
            if !seen_data {
                columns.push(Vec::new());
            } else if k >= columns.len() {
                ogeom_bail!(
                    Construction,
                    "line {number} has more columns than the {} the first point has",
                    columns.len()
                );
            }
            columns[k].push(value);
            width = k + 1;
        }
        if seen_data && width != columns.len() {
            ogeom_bail!(
                Construction,
                "line {number} has {width} columns where the first point has {}",
                columns.len()
            );
        }
        seen_data = true;
    }
    if let Some(left) = owed.filter(|&n| n > 0) {
        ogeom_bail!(
            Construction,
            "the PTS ends with its last scan {left} points short of its count"
        );
    }
    Ok(Table {
        columns,
        header: header_line,
    })
}

/// Whether a column triple reads as 0-to-255 colour channels rather than a
/// normal: every value an integer in `[0, 255]`, and one of them past the
/// `1` no component of a unit normal reaches.
fn colour_like(columns: &[Vec<f64>]) -> bool {
    let channel = |v: &f64| v.fract() == 0.0 && (0.0..=255.0).contains(v);
    columns.iter().all(|c| c.iter().all(channel))
        && columns.iter().any(|c| c.iter().any(|&v| v > 1.0))
}

/// The columns a cloud of `width` unnamed columns is read with.
fn inferred(columns: &[Vec<f64>]) -> OgeomResult<Vec<XyzColumn>> {
    use XyzColumn::{Blue, Green, NormalX, NormalY, NormalZ, Red, X, Y, Z};
    let intensity = || XyzColumn::Scalar("intensity".to_owned());
    let colour = [Red, Green, Blue];
    let normal = [NormalX, NormalY, NormalZ];
    let mut layout = vec![X, Y, Z];
    match columns.len() {
        3 => {}
        4 => layout.push(intensity()),
        6 => layout.extend(if colour_like(&columns[3..6]) {
            colour
        } else {
            normal
        }),
        7 => {
            layout.push(intensity());
            layout.extend(colour);
        }
        9 => match (colour_like(&columns[3..6]), colour_like(&columns[6..9])) {
            (true, false) => {
                layout.extend(colour);
                layout.extend(normal);
            }
            (false, true) => {
                layout.extend(normal);
                layout.extend(colour);
            }
            _ => ogeom_bail!(
                Construction,
                "a cloud of nine columns reads as colour and normals only when one triple is colour; name its columns"
            ),
        },
        n => ogeom_bail!(
            Construction,
            "a cloud of {n} columns does not say what they are; name its columns"
        ),
    }
    Ok(layout)
}

fn assemble(table: Table, named: Option<&[XyzColumn]>) -> OgeomResult<PointCloud> {
    let Table { columns, header } = table;
    if columns.is_empty() {
        ogeom_bail!(Construction, "the file carries no points");
    }
    let layout = match (named, header) {
        (Some(named), _) => named.to_vec(),
        (None, Some(header)) => header,
        (None, None) => inferred(&columns)?,
    };
    if layout.len() != columns.len() {
        ogeom_bail!(
            Construction,
            "the rows have {} columns where {} are named",
            columns.len(),
            layout.len()
        );
    }
    let find = |wanted: &XyzColumn| -> OgeomResult<Option<usize>> {
        let mut at = layout.iter().enumerate().filter(|(_, c)| *c == wanted);
        match (at.next(), at.next()) {
            (Some((i, _)), None) => Ok(Some(i)),
            (None, _) => Ok(None),
            (Some(_), Some(_)) => {
                ogeom_bail!(Construction, "the column {} is named twice", wanted.word())
            }
        }
    };
    let triple = |names: [XyzColumn; 3], what: &str| -> OgeomResult<Option<[usize; 3]>> {
        match (find(&names[0])?, find(&names[1])?, find(&names[2])?) {
            (Some(a), Some(b), Some(c)) => Ok(Some([a, b, c])),
            (None, None, None) => Ok(None),
            _ => ogeom_bail!(Construction, "the columns name part of a {what}"),
        }
    };
    let Some(xyz) = triple([XyzColumn::X, XyzColumn::Y, XyzColumn::Z], "point")? else {
        ogeom_bail!(Construction, "no columns are named x, y and z");
    };
    let nxyz = triple(
        [XyzColumn::NormalX, XyzColumn::NormalY, XyzColumn::NormalZ],
        "normal",
    )?;
    let rgb = triple(
        [XyzColumn::Red, XyzColumn::Green, XyzColumn::Blue],
        "colour",
    )?;
    let n = columns[0].len();
    let points = (0..n)
        .map(|i| Point::new(columns[xyz[0]][i], columns[xyz[1]][i], columns[xyz[2]][i]))
        .collect();
    let normals = nxyz.map(|[a, b, c]| {
        (0..n)
            .map(|i| Vector::new(columns[a][i], columns[b][i], columns[c][i]))
            .collect()
    });
    let colours = match rgb {
        Some(channels) => {
            for &c in &channels {
                if let Some(v) = columns[c].iter().find(|v| !(0.0..=255.0).contains(*v)) {
                    ogeom_bail!(Construction, "a colour channel reads {v}, outside 0 to 255");
                }
            }
            Some(
                (0..n)
                    .map(|i| {
                        let [r, g, b] = channels.map(|c| columns[c][i] / 255.0);
                        [r, g, b, 1.0]
                    })
                    .collect(),
            )
        }
        None => None,
    };
    let mut scalars: Vec<Scalar> = Vec::new();
    let mut columns = columns;
    for (k, column) in layout.iter().enumerate() {
        if let XyzColumn::Scalar(name) = column {
            if scalars.iter().any(|s| &s.name == name) {
                ogeom_bail!(Construction, "the column {name} is named twice");
            }
            scalars.push(Scalar {
                name: name.clone(),
                values: std::mem::take(&mut columns[k]),
            });
        }
    }
    Ok(PointCloud {
        points,
        normals,
        colours,
        scalars,
    })
}

/// Read an XYZ text cloud.
///
/// A row per point; blank lines and lines opening with `#` or `//` are
/// passed over, but for the first such line before the data when it names
/// the columns, `x y z` first: `nx ny nz` are a normal, `r g b` (or `red
/// green blue`) a colour from 0 to 255, and any other word a scalar of
/// that name. `columns`, when given, names them instead. Unnamed, the
/// column count decides:
///
/// | columns | read as |
/// |---------|---------|
/// | 3 | `x y z` |
/// | 4 | `x y z`, then a scalar `intensity` |
/// | 6 | `x y z`, then a colour if every value of the last three is an integer from 0 to 255 and one is past 1, a normal otherwise |
/// | 7 | `x y z`, `intensity`, then a colour, as PTS lays them out |
/// | 9 | `x y z`, then a colour and a normal, the colour being whichever triple reads as one by the rule above |
///
/// Every other count is refused, as is a nine-column file where neither
/// triple, or both, reads as colour.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) for a
/// field that is not a finite number, rows of differing widths, a colour
/// channel outside 0 to 255, columns naming part of a point, normal or
/// colour, or a layout the rule above leaves open.
pub fn read_xyz(text: &str, columns: Option<&[XyzColumn]>) -> OgeomResult<PointCloud> {
    assemble(table(text, false)?, columns)
}

/// Read a PTS text cloud: XYZ rows, each scan opened by a line holding its
/// point count, the rows commonly `x y z intensity red green blue`.
///
/// The columns are read as [`read_xyz`] reads them; several scans in one
/// file are read one after the other into one cloud.
///
/// # Errors
///
/// As [`read_xyz`], and for a point before any count, a scan with more or
/// fewer rows than its count, or a count line that is not a count.
pub fn read_pts(text: &str, columns: Option<&[XyzColumn]>) -> OgeomResult<PointCloud> {
    assemble(table(text, true)?, columns)
}

/// Write a cloud as XYZ text.
///
/// A header line naming the columns, then a row per point: `x y z`, the
/// normal, the colour as 0-to-255 integers, then each scalar, every number
/// in the shortest form that reads back as the same `f64`. Alpha is not
/// written; the format has no column for it. [`read_xyz`] reads the
/// header, so the file comes back as the cloud it was, colours rounded to
/// the nearest 255th.
///
/// # Errors
///
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) if an
/// attribute does not have one entry per point (see [`PointCloud::check`]),
/// and [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for a scalar
/// name that is empty, holds a separator or a `#`, or is one the header
/// gives to a coordinate, normal or colour.
pub fn write_xyz(cloud: &PointCloud) -> OgeomResult<String> {
    cloud.check()?;
    let mut out = String::from("# x y z");
    if cloud.normals.is_some() {
        out.push_str(" nx ny nz");
    }
    if cloud.colours.is_some() {
        out.push_str(" red green blue");
    }
    for s in &cloud.scalars {
        let separated = s
            .name
            .chars()
            .any(|c| c.is_whitespace() || c == ',' || c == ';' || c == '#');
        if s.name.is_empty()
            || separated
            || !matches!(XyzColumn::named(&s.name), XyzColumn::Scalar(_))
        {
            ogeom_bail!(
                Domain,
                "the scalar name {:?} cannot head an XYZ column",
                s.name
            );
        }
        let _ = write!(out, " {}", s.name);
    }
    out.push('\n');
    for (i, p) in cloud.points.iter().enumerate() {
        let _ = write!(out, "{} {} {}", p.x, p.y, p.z);
        if let Some(normals) = &cloud.normals {
            let n = normals[i];
            let _ = write!(out, " {} {} {}", n.x, n.y, n.z);
        }
        if let Some(colours) = &cloud.colours {
            let [r, g, b] = [0, 1, 2].map(|k| channel_byte(colours[i][k]));
            let _ = write!(out, " {r} {g} {b}");
        }
        for s in &cloud.scalars {
            let _ = write!(out, " {}", s.values[i]);
        }
        out.push('\n');
    }
    Ok(out)
}

/// A colour channel in `[0, 1]` as the nearest of 0 to 255.
pub(crate) fn channel_byte(value: f64) -> u8 {
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "clamped into 0..=255 first"
    )]
    let byte = (value.clamp(0.0, 1.0) * 255.0).round() as u8;
    byte
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;

    fn sample_cloud() -> PointCloud {
        PointCloud {
            points: vec![Point::new(0.1, -2.5, 1e-7), Point::new(3.0, 4.0, 5.0)],
            normals: Some(vec![Vector::Z, Vector::new(0.6, 0.0, -0.8)]),
            colours: Some(vec![[1.0, 0.0, 51.0 / 255.0, 1.0], [0.0, 1.0, 0.0, 1.0]]),
            scalars: vec![Scalar {
                name: "deviation".into(),
                values: vec![0.012_345_678_9, -0.25],
            }],
        }
    }

    #[test]
    fn an_xyz_cloud_round_trips_with_its_header() {
        let cloud = sample_cloud();
        let text = write_xyz(&cloud).unwrap();
        assert!(text.starts_with("# x y z nx ny nz red green blue deviation\n"));
        assert_eq!(read_xyz(&text, None).unwrap(), cloud);
    }

    #[test]
    fn unnamed_columns_are_read_by_their_count() {
        let rgb = read_xyz("1,2,3,255,0,10\n4,5,6,0,128,0\n", None).unwrap();
        assert_eq!(rgb.colours.unwrap()[0], [1.0, 0.0, 10.0 / 255.0, 1.0]);
        assert!(rgb.normals.is_none());
        let normal = read_xyz("1 2 3 0 0 1\n4 5 6 0 1 0\n", None).unwrap();
        assert_eq!(normal.normals.unwrap()[1], Vector::Y);
        let intensity = read_xyz("1 2 3 -40\n", None).unwrap();
        assert_eq!(intensity.scalar("intensity"), Some(&[-40.0][..]));
        let nine = read_xyz("0 0 0 0.6 0 0.8 10 20 30\n", None).unwrap();
        assert_eq!(nine.normals.unwrap()[0], Vector::new(0.6, 0.0, 0.8));
        assert_eq!(nine.colours.unwrap()[0][1], 20.0 / 255.0);
        assert!(read_xyz("1 2 3 4 5\n", None).is_err());
        assert!(read_xyz("0 0 0 10 20 30 10 20 30\n", None).is_err());
        let named = read_xyz(
            "1 2 3 9 7\n",
            Some(&[
                XyzColumn::Skip,
                XyzColumn::X,
                XyzColumn::Y,
                XyzColumn::Z,
                XyzColumn::Scalar("quality".into()),
            ]),
        )
        .unwrap();
        assert_eq!(named.points[0], Point::new(2.0, 3.0, 9.0));
        assert_eq!(named.scalar("quality"), Some(&[7.0][..]));
    }

    #[test]
    fn a_pts_reads_every_scan_and_holds_each_to_its_count() {
        let text = "2\n0 0 0 -100 255 0 0\n1 0 0 200 0 255 0\n1\n0 1 0 0 0 0 255\n";
        let cloud = read_pts(text, None).unwrap();
        assert_eq!(cloud.len(), 3);
        assert_eq!(cloud.scalar("intensity"), Some(&[-100.0, 200.0, 0.0][..]));
        assert_eq!(cloud.colours.unwrap()[2], [0.0, 0.0, 1.0, 1.0]);
        assert!(read_pts("3\n0 0 0\n1 1 1\n", None).is_err());
        assert!(read_pts("1\n0 0 0\n1 1 1\n", None).is_err());
        assert!(read_pts("0 0 0\n", None).is_err());
    }

    #[test]
    fn garbage_and_truncation_are_refused_by_name() {
        for bad in [
            "",
            "# x y z\n",
            "1 2 nan\n",
            "1 2 3\n4 5\n",
            "1 2 3\n4 5 6 7\n",
            "1 2 x\n",
            "1 2 3 0 300 0 1\n",
            "\u{0}\u{1}\u{2}",
        ] {
            assert!(read_xyz(bad, None).is_err(), "{bad:?}");
        }
        let text = write_xyz(&sample_cloud()).unwrap();
        for cut in 1..text.len() {
            // Every prefix either reads or is refused; none panics.
            let _ = read_xyz(&text[..cut], None);
        }
        let mut bad_name = sample_cloud();
        bad_name.scalars[0].name = "red".into();
        assert!(write_xyz(&bad_name).is_err());
        let mut short = sample_cloud();
        short.normals = Some(vec![Vector::Z]);
        assert!(write_xyz(&short).is_err());
    }
}
