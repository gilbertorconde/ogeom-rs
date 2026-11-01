//! Texture images, decoded for the colours a mesh samples from them.
//!
//! PNG is decoded here, over the inflater the 3MF reader keeps: chunks
//! checked against their checksums, every colour type and bit depth the
//! format defines, the five scanline filters, Adam7 interlacing, and
//! transparency from `tRNS`. JPEG goes through `zune-jpeg`. Either comes
//! back as 8-bit RGBA, rows from the top.
//!
//! A mesh is sampled at its texture coordinates by the nearest texel, the
//! one whose square holds the coordinate, so a texture painted in flat
//! regions gives back exactly the colours it was painted in.

use ogeom_core::{OgeomResult, ogeom_bail};

/// The longest side a texture may have: 16 384 texels, the largest a
/// graphics device is required to accept. A header promising more is
/// refused before anything is allocated for it.
pub(crate) const MOST_SIDE: usize = 1 << 14;

/// A decoded image: `width × height` texels, row by row from the top.
#[derive(Debug, Clone)]
pub(crate) struct Image {
    pub(crate) width: usize,
    pub(crate) height: usize,
    pub(crate) texels: Vec<[u8; 4]>,
}

/// What a texture coordinate outside `[0, 1]` reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Wrap {
    /// The texture tiles.
    Repeat,
    /// The texture tiles, every other copy mirrored.
    Mirror,
    /// The edge texel carries on.
    Clamp,
}

impl Image {
    /// The colour at `(u, v)`, RGBA in `[0, 1]`, where `v = 0` is the
    /// bottom row, as OBJ and 3MF both place it. A coordinate within
    /// `[0, 1]` is read where it lies, its ends included, whatever the wrap.
    pub(crate) fn sample(&self, u: f64, v: f64, wrap: [Wrap; 2]) -> [f64; 4] {
        let x = texel(u, self.width, wrap[0]);
        let y = self.height - 1 - texel(v, self.height, wrap[1]);
        self.texels[y * self.width + x].map(|c| f64::from(c) / 255.0)
    }
}

/// The texel along one axis of `n` that holds `t`.
fn texel(t: f64, n: usize, wrap: Wrap) -> usize {
    let t = if (0.0..=1.0).contains(&t) {
        t
    } else {
        match wrap {
            Wrap::Repeat => t - t.floor(),
            Wrap::Mirror => {
                let m = t.rem_euclid(2.0);
                if m > 1.0 { 2.0 - m } else { m }
            }
            Wrap::Clamp => t.clamp(0.0, 1.0),
        }
    };
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        clippy::cast_precision_loss,
        reason = "t is in [0, 1] and n at most MOST_SIDE; a cast saturates"
    )]
    let at = (t * n as f64).floor() as usize;
    at.min(n - 1)
}

/// Decode a PNG or a JPEG, told apart by their signatures.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) for
/// bytes that are neither, a file that is damaged or cut short, or an image
/// wider or taller than [`MOST_SIDE`].
pub(crate) fn decode(bytes: &[u8]) -> OgeomResult<Image> {
    if bytes.starts_with(b"\x89PNG\r\n\x1a\n") {
        png(bytes)
    } else if bytes.starts_with(&[0xFF, 0xD8]) {
        jpeg(bytes)
    } else {
        ogeom_bail!(Construction, "the image is neither a PNG nor a JPEG")
    }
}

/// The texels the textures of one file may hold together: as many as one
/// texture of the longest side. A file whose textures would hold more has
/// the rest refused by name, each before it is decoded.
pub(crate) const MOST_TEXELS: usize = MOST_SIDE * MOST_SIDE;

/// What is left of one file's texel allowance.
pub(crate) struct Budget(usize);

impl Budget {
    pub(crate) fn new() -> Self {
        Self(MOST_TEXELS)
    }

    /// [`decode`], refused before decoding if the image would take the
    /// file's textures past [`MOST_TEXELS`].
    pub(crate) fn decode(&mut self, bytes: &[u8]) -> OgeomResult<Image> {
        let (width, height) = dimensions(bytes)?;
        check_size(width, height)?;
        let texels = width * height;
        if texels > self.0 {
            ogeom_bail!(
                Construction,
                "the image of {width} by {height} takes the file's textures past {MOST_TEXELS} texels"
            );
        }
        let image = decode(bytes)?;
        self.0 -= texels;
        Ok(image)
    }
}

/// The width and height an image's header declares.
fn dimensions(bytes: &[u8]) -> OgeomResult<(usize, usize)> {
    if bytes.starts_with(b"\x89PNG\r\n\x1a\n") {
        if bytes.get(12..16) != Some(b"IHDR") {
            ogeom_bail!(Construction, "the PNG does not open with its header chunk");
        }
        return Ok((be32(bytes, 16)? as usize, be32(bytes, 20)? as usize));
    }
    if bytes.starts_with(&[0xFF, 0xD8]) {
        let mut decoder =
            zune_jpeg::JpegDecoder::new(zune_jpeg::zune_core::bytestream::ZCursor::new(bytes));
        if let Err(e) = decoder.decode_headers() {
            ogeom_bail!(Construction, "the JPEG cannot be read: {e}");
        }
        let Some(info) = decoder.info() else {
            ogeom_bail!(Construction, "the JPEG has no frame header");
        };
        return Ok((usize::from(info.width), usize::from(info.height)));
    }
    ogeom_bail!(Construction, "the image is neither a PNG nor a JPEG")
}

fn check_size(width: usize, height: usize) -> OgeomResult<()> {
    if width == 0 || height == 0 {
        ogeom_bail!(
            Construction,
            "the image is {width} by {height}, which holds nothing"
        );
    }
    if width > MOST_SIDE || height > MOST_SIDE {
        ogeom_bail!(
            Construction,
            "the image is {width} by {height}, past the {MOST_SIDE} a side a texture may have"
        );
    }
    Ok(())
}

fn jpeg(bytes: &[u8]) -> OgeomResult<Image> {
    use zune_jpeg::zune_core::bytestream::ZCursor;
    use zune_jpeg::zune_core::colorspace::ColorSpace;
    use zune_jpeg::zune_core::options::DecoderOptions;
    let options = DecoderOptions::default()
        .jpeg_set_out_colorspace(ColorSpace::RGBA)
        .set_max_width(MOST_SIDE)
        .set_max_height(MOST_SIDE);
    let mut decoder = zune_jpeg::JpegDecoder::new_with_options(ZCursor::new(bytes), options);
    if let Err(e) = decoder.decode_headers() {
        ogeom_bail!(Construction, "the JPEG cannot be read: {e}");
    }
    let Some(info) = decoder.info() else {
        ogeom_bail!(Construction, "the JPEG has no frame header");
    };
    let (width, height) = (usize::from(info.width), usize::from(info.height));
    check_size(width, height)?;
    let pixels = match decoder.decode() {
        Ok(p) => p,
        Err(e) => ogeom_bail!(Construction, "the JPEG cannot be read: {e}"),
    };
    if pixels.len() != width * height * 4 {
        ogeom_bail!(
            Construction,
            "the JPEG decodes to {} bytes for {width} by {height} texels",
            pixels.len()
        );
    }
    let texels = pixels.as_chunks::<4>().0.to_vec();
    Ok(Image {
        width,
        height,
        texels,
    })
}

// --- PNG ---------------------------------------------------------------------

/// The header's fields that shape the pixel data.
struct Header {
    width: usize,
    height: usize,
    depth: u8,
    colour: u8,
    interlaced: bool,
}

impl Header {
    fn channels(&self) -> usize {
        match self.colour {
            2 => 3,
            4 => 2,
            6 => 4,
            _ => 1,
        }
    }

    fn bits_per_pixel(&self) -> usize {
        self.channels() * usize::from(self.depth)
    }

    /// Bytes in one filtered row of `width` pixels, its filter byte aside.
    fn row_bytes(&self, width: usize) -> usize {
        (width * self.bits_per_pixel()).div_ceil(8)
    }
}

/// The seven Adam7 passes: start column, start row, column step, row step.
const ADAM7: [(usize, usize, usize, usize); 7] = [
    (0, 0, 8, 8),
    (4, 0, 8, 8),
    (0, 4, 4, 8),
    (2, 0, 4, 4),
    (0, 2, 2, 4),
    (1, 0, 2, 2),
    (0, 1, 1, 2),
];

/// Each pass of the image: `(x0, y0, dx, dy, columns, rows)`; one pass
/// covering every pixel when the image is not interlaced.
fn passes(header: &Header) -> Vec<(usize, usize, usize, usize, usize, usize)> {
    if !header.interlaced {
        return vec![(0, 0, 1, 1, header.width, header.height)];
    }
    ADAM7
        .iter()
        .map(|&(x0, y0, dx, dy)| {
            let columns = header.width.saturating_sub(x0).div_ceil(dx);
            let rows = header.height.saturating_sub(y0).div_ceil(dy);
            (x0, y0, dx, dy, columns, rows)
        })
        .collect()
}

fn be32(bytes: &[u8], at: usize) -> OgeomResult<u32> {
    match at.checked_add(4).and_then(|end| bytes.get(at..end)) {
        Some(b) => Ok(u32::from_be_bytes([b[0], b[1], b[2], b[3]])),
        None => ogeom_bail!(Construction, "the PNG is cut short"),
    }
}

fn png(bytes: &[u8]) -> OgeomResult<Image> {
    let mut at = 8;
    let mut header: Option<Header> = None;
    let mut palette: Vec<[u8; 4]> = Vec::new();
    let mut transparent: Option<Vec<u8>> = None;
    let mut data: Vec<u8> = Vec::new();
    loop {
        let length = be32(bytes, at)? as usize;
        let Some(chunk) = at
            .checked_add(12 + length)
            .and_then(|end| bytes.get(at + 4..end))
        else {
            ogeom_bail!(Construction, "the PNG is cut short");
        };
        let (kind, rest) = chunk.split_at(4);
        let (body, crc) = rest.split_at(length);
        if crate::threemf::crc32(&chunk[..4 + length]) != be32(crc, 0)? {
            ogeom_bail!(
                Construction,
                "the PNG chunk {} does not match its checksum",
                String::from_utf8_lossy(kind)
            );
        }
        at += 12 + length;
        match kind {
            b"IHDR" => {
                if body.len() != 13 {
                    ogeom_bail!(
                        Construction,
                        "the PNG header is {} bytes, not 13",
                        body.len()
                    );
                }
                let width = be32(body, 0)? as usize;
                let height = be32(body, 4)? as usize;
                check_size(width, height)?;
                let (depth, colour) = (body[8], body[9]);
                let allowed: &[u8] = match colour {
                    0 => &[1, 2, 4, 8, 16],
                    3 => &[1, 2, 4, 8],
                    2 | 4 | 6 => &[8, 16],
                    _ => ogeom_bail!(
                        Construction,
                        "the PNG colour type {colour} is not one the format defines"
                    ),
                };
                if !allowed.contains(&depth) {
                    ogeom_bail!(
                        Construction,
                        "the PNG colour type {colour} has no bit depth {depth}"
                    );
                }
                if body[10] != 0 || body[11] != 0 || body[12] > 1 {
                    ogeom_bail!(
                        Construction,
                        "the PNG names a compression, filter or interlace method the format does not define"
                    );
                }
                header = Some(Header {
                    width,
                    height,
                    depth,
                    colour,
                    interlaced: body[12] == 1,
                });
            }
            b"PLTE" => {
                if body.len() % 3 != 0 || body.len() > 3 * 256 {
                    ogeom_bail!(
                        Construction,
                        "the PNG palette is not whole entries of three bytes"
                    );
                }
                palette = body
                    .as_chunks::<3>()
                    .0
                    .iter()
                    .map(|&[r, g, b]| [r, g, b, 255])
                    .collect();
            }
            b"tRNS" => transparent = Some(body.to_vec()),
            b"IDAT" => data.extend_from_slice(body),
            b"IEND" => break,
            _ if kind[0].is_ascii_uppercase() => ogeom_bail!(
                Construction,
                "the PNG carries the critical chunk {}, which is not read",
                String::from_utf8_lossy(kind)
            ),
            _ => {}
        }
    }
    let Some(header) = header else {
        ogeom_bail!(Construction, "the PNG has no header chunk");
    };
    if header.colour == 3 {
        if palette.is_empty() {
            ogeom_bail!(Construction, "the PNG is indexed and has no palette");
        }
        if let Some(alpha) = &transparent {
            for (entry, &a) in palette.iter_mut().zip(alpha) {
                entry[3] = a;
            }
        }
    }

    // A zlib stream: a two-byte header, then raw DEFLATE.
    if data.len() < 2
        || data[0] & 0x0F != 8
        || (u16::from(data[0]) << 8 | u16::from(data[1])) % 31 != 0
        || data[1] & 0x20 != 0
    {
        ogeom_bail!(Construction, "the PNG's image data is not a zlib stream");
    }
    let layout = passes(&header);
    let expected: usize = layout
        .iter()
        .filter(|p| p.4 > 0 && p.5 > 0)
        .map(|p| p.5 * (1 + header.row_bytes(p.4)))
        .sum();
    let raw = crate::inflate::inflate(&data[2..], expected)?;

    let mut texels = vec![[0_u8; 4]; header.width * header.height];
    let stride = header.bits_per_pixel().div_ceil(8);
    let mut offset = 0;
    for (x0, y0, dx, dy, columns, rows) in layout {
        if columns == 0 || rows == 0 {
            continue;
        }
        let row_bytes = header.row_bytes(columns);
        let mut previous = vec![0_u8; row_bytes];
        let mut current = vec![0_u8; row_bytes];
        for row in 0..rows {
            let filter = raw[offset];
            current.copy_from_slice(&raw[offset + 1..offset + 1 + row_bytes]);
            offset += 1 + row_bytes;
            unfilter(filter, &mut current, &previous, stride)?;
            for column in 0..columns {
                let texel = pixel(&header, &current, column, &palette, transparent.as_deref())?;
                texels[(y0 + row * dy) * header.width + x0 + column * dx] = texel;
            }
            std::mem::swap(&mut previous, &mut current);
        }
    }
    Ok(Image {
        width: header.width,
        height: header.height,
        texels,
    })
}

/// Undo one scanline's filter in place, `previous` the row above it
/// unfiltered, `stride` the bytes per pixel (at least one).
fn unfilter(filter: u8, row: &mut [u8], previous: &[u8], stride: usize) -> OgeomResult<()> {
    for i in 0..row.len() {
        let left = if i >= stride { row[i - stride] } else { 0 };
        let up = previous[i];
        let corner = if i >= stride { previous[i - stride] } else { 0 };
        let predicted = match filter {
            0 => 0,
            1 => left,
            2 => up,
            #[allow(clippy::cast_possible_truncation, reason = "a mean of two bytes")]
            3 => ((u16::from(left) + u16::from(up)) / 2) as u8,
            4 => paeth(left, up, corner),
            _ => ogeom_bail!(Construction, "the PNG has a row with filter type {filter}"),
        };
        row[i] = row[i].wrapping_add(predicted);
    }
    Ok(())
}

fn paeth(a: u8, b: u8, c: u8) -> u8 {
    let p = i16::from(a) + i16::from(b) - i16::from(c);
    let (pa, pb, pc) = (
        (p - i16::from(a)).abs(),
        (p - i16::from(b)).abs(),
        (p - i16::from(c)).abs(),
    );
    if pa <= pb && pa <= pc {
        a
    } else if pb <= pc {
        b
    } else {
        c
    }
}

/// Sample `k` of a row, at the image's bit depth, unscaled.
fn sample(row: &[u8], k: usize, depth: u8) -> u16 {
    match depth {
        16 => u16::from_be_bytes([row[2 * k], row[2 * k + 1]]),
        8 => u16::from(row[k]),
        _ => {
            let bits = usize::from(depth);
            let byte = row[k * bits / 8];
            let shift = 8 - bits - (k * bits) % 8;
            u16::from(byte >> shift) & ((1 << bits) - 1)
        }
    }
}

/// A sample scaled to eight bits, rounded.
fn eight(value: u16, depth: u8) -> u8 {
    let most = (1_u32 << depth) - 1;
    #[allow(clippy::cast_possible_truncation, reason = "at most 255")]
    let scaled = ((u32::from(value) * 255 + most / 2) / most) as u8;
    scaled
}

fn pixel(
    header: &Header,
    row: &[u8],
    column: usize,
    palette: &[[u8; 4]],
    transparent: Option<&[u8]>,
) -> OgeomResult<[u8; 4]> {
    let channels = header.channels();
    let depth = header.depth;
    let at = |c: usize| sample(row, column * channels + c, depth);
    let key = |i: usize| -> Option<u16> {
        let t = transparent?;
        Some(u16::from_be_bytes([*t.get(2 * i)?, *t.get(2 * i + 1)?]))
    };
    Ok(match header.colour {
        0 => {
            let grey = at(0);
            let g = eight(grey, depth);
            [g, g, g, if key(0) == Some(grey) { 0 } else { 255 }]
        }
        2 => {
            let (r, g, b) = (at(0), at(1), at(2));
            let opaque = !(key(0) == Some(r) && key(1) == Some(g) && key(2) == Some(b));
            [
                eight(r, depth),
                eight(g, depth),
                eight(b, depth),
                if opaque { 255 } else { 0 },
            ]
        }
        3 => {
            let index = usize::from(at(0));
            let Some(&entry) = palette.get(index) else {
                ogeom_bail!(
                    Construction,
                    "the PNG names palette entry {index} of {}",
                    palette.len()
                );
            };
            entry
        }
        4 => {
            let g = eight(at(0), depth);
            [g, g, g, eight(at(1), depth)]
        }
        _ => [
            eight(at(0), depth),
            eight(at(1), depth),
            eight(at(2), depth),
            eight(at(3), depth),
        ],
    })
}

/// A PNG of 8-bit RGBA texels, unfiltered and stored uncompressed: what
/// the tests here and in the readers sample from.
#[cfg(test)]
pub(crate) fn encode_png(width: usize, height: usize, texels: &[[u8; 4]]) -> Vec<u8> {
    fn chunk(out: &mut Vec<u8>, kind: &[u8], body: &[u8]) {
        out.extend_from_slice(&u32::try_from(body.len()).unwrap_or(0).to_be_bytes());
        let mut typed = kind.to_vec();
        typed.extend_from_slice(body);
        out.extend_from_slice(&typed);
        out.extend_from_slice(&crate::threemf::crc32(&typed).to_be_bytes());
    }
    let mut raw = Vec::new();
    for row in texels.chunks(width) {
        raw.push(0);
        for t in row {
            raw.extend_from_slice(t);
        }
    }
    // zlib header, then stored DEFLATE blocks of at most 65 535 bytes.
    let mut z = vec![0x78, 0x01];
    let blocks: Vec<&[u8]> = raw.chunks(0xFFFF).collect();
    for (i, block) in blocks.iter().enumerate() {
        z.push(u8::from(i + 1 == blocks.len()));
        let len = u16::try_from(block.len()).unwrap_or(0);
        z.extend_from_slice(&len.to_le_bytes());
        z.extend_from_slice(&(!len).to_le_bytes());
        z.extend_from_slice(block);
    }
    z.extend_from_slice(&[0, 0, 0, 0]);
    let mut out = b"\x89PNG\r\n\x1a\n".to_vec();
    let mut ihdr = Vec::new();
    ihdr.extend_from_slice(&u32::try_from(width).unwrap_or(0).to_be_bytes());
    ihdr.extend_from_slice(&u32::try_from(height).unwrap_or(0).to_be_bytes());
    ihdr.extend_from_slice(&[8, 6, 0, 0, 0]);
    chunk(&mut out, b"IHDR", &ihdr);
    chunk(&mut out, b"IDAT", &z);
    chunk(&mut out, b"IEND", &[]);
    out
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;

    #[test]
    fn a_png_reads_back_texel_for_texel() {
        let texels: Vec<[u8; 4]> = (0..12_u8).map(|i| [i * 20, 255 - i, i, 200]).collect();
        let image = decode(&encode_png(4, 3, &texels)).unwrap();
        assert_eq!((image.width, image.height), (4, 3));
        assert_eq!(image.texels, texels);
        // v = 0 is the bottom row; each corner reads its own texel.
        assert_eq!(
            image.sample(0.0, 1.0, [Wrap::Clamp; 2]),
            texels[0].map(|c| f64::from(c) / 255.0)
        );
        assert_eq!(
            image.sample(1.0, 0.0, [Wrap::Clamp; 2]),
            texels[11].map(|c| f64::from(c) / 255.0)
        );
        assert_eq!(
            image.sample(1.1, -0.9, [Wrap::Repeat; 2]),
            image.sample(0.1, 0.1, [Wrap::Clamp; 2])
        );
        assert_eq!(
            image.sample(1.1, 0.5, [Wrap::Mirror; 2]),
            image.sample(0.9, 0.5, [Wrap::Clamp; 2])
        );
    }

    /// Every filter type and an interlaced, palette-indexed, two-bit image,
    /// written out by hand: the filters are undone against their
    /// neighbours, and the passes land on their own pixels.
    #[test]
    fn filters_palettes_and_interlacing_decode() {
        // A 3 × 1 RGB image, each row of one pixel filtered differently.
        let header = Header {
            width: 1,
            height: 5,
            depth: 8,
            colour: 2,
            interlaced: false,
        };
        let rows: [[u8; 3]; 5] = [
            [10, 20, 30],
            [40, 50, 60],
            [5, 5, 5],
            [200, 100, 0],
            [1, 2, 3],
        ];
        let mut previous = [0_u8; 3];
        for (filter, row) in rows.iter().enumerate() {
            let filter = u8::try_from(filter).unwrap();
            // Filter this row against the one above, then undo it.
            let mut encoded = *row;
            for i in 0..3 {
                let up = previous[i];
                let predicted = match filter {
                    0 | 1 => 0, // one pixel per row: nothing to its left
                    2 => up,
                    3 => up / 2,
                    _ => paeth(0, up, 0),
                };
                encoded[i] = row[i].wrapping_sub(predicted);
            }
            unfilter(filter, &mut encoded, &previous, 3).unwrap();
            assert_eq!(&encoded, row);
            assert_eq!(
                pixel(&header, &encoded, 0, &[], None).unwrap()[..3],
                row[..]
            );
            previous = *row;
        }
        assert!(unfilter(5, &mut [0], &[0], 1).is_err());

        // A 2-bit indexed sample, the third of a byte.
        let palette = [[0, 0, 0, 255], [9, 9, 9, 255], [1, 2, 3, 4], [7, 7, 7, 255]];
        let indexed = Header {
            width: 4,
            height: 1,
            depth: 2,
            colour: 3,
            interlaced: true,
        };
        assert_eq!(
            pixel(&indexed, &[0b00_01_10_11], 2, &palette, None).unwrap(),
            [1, 2, 3, 4]
        );
        // Adam7 covers every pixel of a 5 × 3 image once.
        let mut seen = [[0; 5]; 3];
        for (x0, y0, dx, dy, columns, rows) in passes(&Header {
            width: 5,
            height: 3,
            ..indexed
        }) {
            for r in 0..rows {
                for c in 0..columns {
                    seen[y0 + r * dy][x0 + c * dx] += 1;
                }
            }
        }
        assert!(seen.iter().flatten().all(|&n| n == 1));
    }

    #[test]
    fn damaged_and_outsized_images_are_refused_by_name() {
        let png = encode_png(2, 2, &[[1, 2, 3, 4]; 4]);
        for cut in [9, 20, 33, png.len() - 1] {
            assert!(decode(&png[..cut]).is_err(), "cut at {cut}");
        }
        let mut damaged = png.clone();
        damaged[45] ^= 0x40;
        assert!(decode(&damaged).is_err());
        let huge = encode_png(1, 1, &[[0; 4]]);
        let mut huge = huge.clone();
        // Widen the header to 2^20 texels and mend its checksum.
        huge[16..20].copy_from_slice(&(1_u32 << 20).to_be_bytes());
        let crc = crate::threemf::crc32(&huge[12..29]);
        huge[29..33].copy_from_slice(&crc.to_be_bytes());
        let error = decode(&huge).unwrap_err().to_string();
        assert!(error.contains("a side"), "{error}");
        // A file's textures share one allowance, checked from the header
        // before anything is decoded.
        let mut budget = Budget(5);
        assert!(budget.decode(&png).is_ok());
        let error = budget.decode(&png).unwrap_err().to_string();
        assert!(error.contains("past 268435456 texels"), "{error}");
        assert!(decode(b"GIF89a").is_err());
        assert!(decode(&[0xFF, 0xD8, 0xFF, 0xE0, 0, 16]).is_err());
    }
}
