//! The ISO 10303-21 exchange structure, parsed but not yet interpreted.
//!
//! Part 21 is a syntax, not a schema: `#12 = CIRCLE('', #10, 5.0);` says an
//! instance exists with a keyword and arguments, and what a `CIRCLE` *means*
//! is the reader's business, not the parser's. This module turns the text
//! into a map from instance number to typed argument trees and nothing more —
//! which is what lets the reader say precisely which entities it understood
//! and which it deliberately walked past.

use ogeom_core::{OgeomResult, ogeom_bail};
use std::collections::HashMap;

/// One argument of an entity instance.
#[derive(Debug, Clone, PartialEq)]
pub enum Arg {
    /// `$` — no value.
    Null,
    /// `*` — value derivable from the schema, not stated.
    Derived,
    /// `#n` — a reference to another instance.
    Ref(u64),
    /// An integer literal.
    Int(i64),
    /// A real literal.
    Real(f64),
    /// A string literal, with Part 21's quote doubling undone.
    Str(String),
    /// `.NAME.` — an enumeration value, without its dots.
    Enum(String),
    /// A parenthesised list.
    List(Vec<Arg>),
    /// `KEYWORD(...)` in argument position — a typed (select) value.
    Typed(String, Vec<Arg>),
}

impl Arg {
    /// The reference this argument carries, if it is one.
    pub fn reference(&self) -> Option<u64> {
        match self {
            Self::Ref(n) => Some(*n),
            _ => None,
        }
    }

    /// The number this argument carries, integer or real.
    pub fn number(&self) -> Option<f64> {
        match self {
            Self::Int(n) =>
            {
                #[allow(clippy::cast_precision_loss)]
                Some(*n as f64)
            }
            Self::Real(x) => Some(*x),
            _ => None,
        }
    }

    /// The list this argument carries, if it is one.
    pub fn list(&self) -> Option<&[Arg]> {
        match self {
            Self::List(items) => Some(items),
            _ => None,
        }
    }

    /// Whether this is the enumeration value `name`.
    pub fn is_enum(&self, name: &str) -> bool {
        matches!(self, Self::Enum(e) if e == name)
    }
}

/// One instance: usually one keyword with arguments, several for a complex
/// (multi-leaf) instance like `#1 = (A(...) B(...));`.
#[derive(Debug, Clone)]
pub struct Instance {
    /// The parts, in file order.
    pub parts: Vec<(String, Vec<Arg>)>,
}

impl Instance {
    /// The arguments of the part with this keyword, if present.
    pub fn part(&self, keyword: &str) -> Option<&[Arg]> {
        self.parts
            .iter()
            .find(|(k, _)| k == keyword)
            .map(|(_, a)| a.as_slice())
    }

    /// The single keyword of a simple instance.
    pub fn keyword(&self) -> &str {
        self.parts.first().map_or("", |(k, _)| k.as_str())
    }
}

/// A parsed exchange file: the data section as a graph, the header kept as
/// raw instances for whoever wants the schema name or the file's own record
/// of itself.
#[derive(Debug)]
pub struct Exchange {
    /// Header entries, in order.
    pub header: Vec<(String, Vec<Arg>)>,
    /// The data section, by instance number.
    pub data: HashMap<u64, Instance>,
}

/// Parse a Part 21 exchange file.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) on malformed
/// syntax, with the byte offset where reading stopped making sense.
pub fn parse(text: &str) -> OgeomResult<Exchange> {
    let mut p = Parser {
        bytes: text.as_bytes(),
        at: 0,
    };
    p.skip_noise();
    p.expect_keyword("ISO-10303-21")?;
    p.expect(b';')?;

    p.expect_keyword("HEADER")?;
    p.expect(b';')?;
    let mut header = Vec::new();
    loop {
        p.skip_noise();
        if p.peek_keyword("ENDSEC") {
            p.expect_keyword("ENDSEC")?;
            p.expect(b';')?;
            break;
        }
        let keyword = p.keyword()?;
        let args = p.arguments()?;
        p.expect(b';')?;
        header.push((keyword, args));
    }

    p.expect_keyword("DATA")?;
    p.expect(b';')?;
    let mut data = HashMap::new();
    loop {
        p.skip_noise();
        if p.peek_keyword("ENDSEC") {
            p.expect_keyword("ENDSEC")?;
            p.expect(b';')?;
            break;
        }
        p.expect(b'#')?;
        let id = p.integer()?;
        p.expect(b'=')?;
        p.skip_noise();
        let parts = if p.peek(b'(') {
            // A complex instance: parenthesised sequence of parts.
            p.expect(b'(')?;
            let mut parts = Vec::new();
            loop {
                p.skip_noise();
                if p.peek(b')') {
                    p.expect(b')')?;
                    break;
                }
                let keyword = p.keyword()?;
                let args = p.arguments()?;
                parts.push((keyword, args));
            }
            parts
        } else {
            let keyword = p.keyword()?;
            let args = p.arguments()?;
            vec![(keyword, args)]
        };
        p.expect(b';')?;
        #[allow(clippy::cast_sign_loss)]
        data.insert(id as u64, Instance { parts });
    }

    p.expect_keyword("END-ISO-10303-21")?;
    Ok(Exchange { header, data })
}

struct Parser<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Parser<'_> {
    fn skip_noise(&mut self) {
        loop {
            while self.at < self.bytes.len() && self.bytes[self.at].is_ascii_whitespace() {
                self.at += 1;
            }
            if self.at + 1 < self.bytes.len() && &self.bytes[self.at..self.at + 2] == b"/*" {
                self.at += 2;
                while self.at + 1 < self.bytes.len() && &self.bytes[self.at..self.at + 2] != b"*/" {
                    self.at += 1;
                }
                self.at = (self.at + 2).min(self.bytes.len());
                continue;
            }
            break;
        }
    }

    fn peek(&mut self, byte: u8) -> bool {
        self.skip_noise();
        self.bytes.get(self.at) == Some(&byte)
    }

    fn expect(&mut self, byte: u8) -> OgeomResult<()> {
        self.skip_noise();
        if self.bytes.get(self.at) == Some(&byte) {
            self.at += 1;
            return Ok(());
        }
        ogeom_bail!(
            Construction,
            "expected '{}' at byte {} of the exchange file",
            char::from(byte),
            self.at
        );
    }

    fn peek_keyword(&mut self, word: &str) -> bool {
        self.skip_noise();
        let end = self.at + word.len();
        end <= self.bytes.len()
            && &self.bytes[self.at..end] == word.as_bytes()
            && self
                .bytes
                .get(end)
                .is_none_or(|b| !b.is_ascii_alphanumeric() && *b != b'_' && *b != b'-')
    }

    fn expect_keyword(&mut self, word: &str) -> OgeomResult<()> {
        if self.peek_keyword(word) {
            self.at += word.len();
            return Ok(());
        }
        ogeom_bail!(
            Construction,
            "expected '{word}' at byte {} of the exchange file",
            self.at
        );
    }

    fn keyword(&mut self) -> OgeomResult<String> {
        self.skip_noise();
        let start = self.at;
        while self
            .bytes
            .get(self.at)
            .is_some_and(|b| b.is_ascii_alphanumeric() || *b == b'_')
        {
            self.at += 1;
        }
        if self.at == start {
            ogeom_bail!(
                Construction,
                "expected a keyword at byte {} of the exchange file",
                start
            );
        }
        Ok(String::from_utf8_lossy(&self.bytes[start..self.at]).into_owned())
    }

    fn integer(&mut self) -> OgeomResult<i64> {
        self.skip_noise();
        let start = self.at;
        if self.bytes.get(self.at) == Some(&b'-') || self.bytes.get(self.at) == Some(&b'+') {
            self.at += 1;
        }
        while self.bytes.get(self.at).is_some_and(u8::is_ascii_digit) {
            self.at += 1;
        }
        let text = std::str::from_utf8(&self.bytes[start..self.at]).unwrap_or("");
        text.parse().map_err(|_| {
            ogeom_core::ogeom_err!(
                Construction,
                "expected an integer at byte {start} of the exchange file"
            )
        })
    }

    fn arguments(&mut self) -> OgeomResult<Vec<Arg>> {
        self.expect(b'(')?;
        let mut out = Vec::new();
        loop {
            self.skip_noise();
            if self.peek(b')') {
                self.expect(b')')?;
                break;
            }
            out.push(self.argument()?);
            self.skip_noise();
            if self.peek(b',') {
                self.expect(b',')?;
            }
        }
        Ok(out)
    }

    fn argument(&mut self) -> OgeomResult<Arg> {
        self.skip_noise();
        let Some(&byte) = self.bytes.get(self.at) else {
            ogeom_bail!(Construction, "the exchange file ends inside an argument");
        };
        match byte {
            b'$' => {
                self.at += 1;
                Ok(Arg::Null)
            }
            b'*' => {
                self.at += 1;
                Ok(Arg::Derived)
            }
            b'#' => {
                self.at += 1;
                let id = self.integer()?;
                #[allow(clippy::cast_sign_loss)]
                Ok(Arg::Ref(id as u64))
            }
            b'(' => Ok(Arg::List(self.arguments()?)),
            b'\'' => self.string(),
            b'.' => {
                self.at += 1;
                let word = self.keyword()?;
                self.expect(b'.')?;
                Ok(Arg::Enum(word))
            }
            b'-' | b'+' | b'0'..=b'9' => self.number(),
            _ if byte.is_ascii_alphabetic() || byte == b'_' => {
                let keyword = self.keyword()?;
                let args = self.arguments()?;
                Ok(Arg::Typed(keyword, args))
            }
            _ => ogeom_bail!(
                Construction,
                "unexpected '{}' at byte {} of the exchange file",
                char::from(byte),
                self.at
            ),
        }
    }

    fn string(&mut self) -> OgeomResult<Arg> {
        self.expect(b'\'')?;
        let mut out = String::new();
        loop {
            match self.bytes.get(self.at) {
                None => ogeom_bail!(Construction, "the exchange file ends inside a string"),
                Some(b'\'') => {
                    if self.bytes.get(self.at + 1) == Some(&b'\'') {
                        out.push('\'');
                        self.at += 2;
                    } else {
                        self.at += 1;
                        break;
                    }
                }
                Some(&b) => {
                    out.push(char::from(b));
                    self.at += 1;
                }
            }
        }
        Ok(Arg::Str(out))
    }

    fn number(&mut self) -> OgeomResult<Arg> {
        let start = self.at;
        if matches!(self.bytes.get(self.at), Some(b'-' | b'+')) {
            self.at += 1;
        }
        let mut real = false;
        while let Some(&b) = self.bytes.get(self.at) {
            match b {
                b'0'..=b'9' => self.at += 1,
                b'.' => {
                    // A dot starts a real — unless it starts an enumeration
                    // hard against the number, which no real file does.
                    real = true;
                    self.at += 1;
                }
                b'E' | b'e' => {
                    real = true;
                    self.at += 1;
                    if matches!(self.bytes.get(self.at), Some(b'-' | b'+')) {
                        self.at += 1;
                    }
                }
                _ => break,
            }
        }
        let text = std::str::from_utf8(&self.bytes[start..self.at]).unwrap_or("");
        if real {
            text.parse().map(Arg::Real).map_err(|_| {
                ogeom_core::ogeom_err!(
                    Construction,
                    "unreadable real at byte {start} of the exchange file"
                )
            })
        } else {
            text.parse().map(Arg::Int).map_err(|_| {
                ogeom_core::ogeom_err!(
                    Construction,
                    "unreadable integer at byte {start} of the exchange file"
                )
            })
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    const SMALL: &str = "ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('a part'),'2;1');
FILE_NAME('p.stp','2020-01-01',(''),(''),'','','');
FILE_SCHEMA(('AP203'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,1.5,-2.E-3));
#2=DIRECTION('',(0.,0.,1.));
#3=AXIS2_PLACEMENT_3D('',#1,#2,$);
#4=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#5)));
#5=SI_UNIT(.MILLI.,.METRE.);
ENDSEC;
END-ISO-10303-21;
";

    #[test]
    fn a_small_file_parses_into_its_instances() {
        let file = parse(SMALL).unwrap();
        assert_eq!(file.header.len(), 3);
        assert_eq!(file.data.len(), 5);

        let point = &file.data[&1];
        assert_eq!(point.keyword(), "CARTESIAN_POINT");
        let coords = point.parts[0].1[1].list().unwrap();
        assert_eq!(coords[0].number(), Some(0.0));
        assert_eq!(coords[1].number(), Some(1.5));
        assert_eq!(coords[2].number(), Some(-2e-3));

        let placement = &file.data[&3];
        assert_eq!(placement.parts[0].1[1].reference(), Some(1));
        assert_eq!(placement.parts[0].1[3], Arg::Null);

        // The complex instance keeps both parts, each with its own arguments.
        let context = &file.data[&4];
        assert_eq!(context.parts.len(), 2);
        assert!(context.part("GLOBAL_UNIT_ASSIGNED_CONTEXT").is_some());

        let unit = &file.data[&5];
        assert!(unit.parts[0].1[0].is_enum("MILLI"));
    }

    #[test]
    fn strings_undouble_their_quotes() {
        let file = parse("ISO-10303-21;HEADER;ENDSEC;DATA;#1=X('it''s');ENDSEC;END-ISO-10303-21;")
            .unwrap();
        assert_eq!(file.data[&1].parts[0].1[0], Arg::Str("it's".into()));
    }

    #[test]
    fn malformed_files_are_refused_with_a_place() {
        assert!(parse("ISO-10303-21;HEADER;DATA;").is_err());
        assert!(parse("not a step file").is_err());
    }
}
