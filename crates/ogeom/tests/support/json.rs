//! Just enough JSON to read a test's data file: objects, arrays, numbers,
//! strings, `true`, `false` and `null`.

/// A JSON value.
#[derive(Debug, Clone)]
pub enum Json {
    Null,
    Bool,
    Number(f64),
    Text(String),
    List(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl Json {
    /// The value read from `text`; `None` where it is not JSON.
    pub fn parse(text: &str) -> Option<Self> {
        let mut reader = Reader {
            bytes: text.as_bytes(),
            at: 0,
        };
        let value = reader.value()?;
        reader.space();
        (reader.at == reader.bytes.len()).then_some(value)
    }

    /// The member `key` of an object.
    pub fn get(&self, key: &str) -> Option<&Self> {
        match self {
            Self::Object(members) => members.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    pub fn number(&self) -> Option<f64> {
        match self {
            Self::Number(x) => Some(*x),
            _ => None,
        }
    }

    pub fn list(&self) -> &[Self] {
        match self {
            Self::List(items) => items,
            _ => &[],
        }
    }

    /// The numbers of a list, `N` of them.
    pub fn numbers<const N: usize>(&self) -> Option<[f64; N]> {
        let items = self.list();
        if items.len() != N {
            return None;
        }
        let mut out = [0.0; N];
        for (slot, item) in out.iter_mut().zip(items) {
            *slot = item.number()?;
        }
        Some(out)
    }

    /// An object's only member, for an enum written `{"variant": {...}}`.
    pub fn variant(&self) -> Option<(&str, &Self)> {
        match self {
            Self::Object(members) if members.len() == 1 => {
                Some((members[0].0.as_str(), &members[0].1))
            }
            _ => None,
        }
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Reader<'_> {
    fn space(&mut self) {
        while self.bytes.get(self.at).is_some_and(u8::is_ascii_whitespace) {
            self.at += 1;
        }
    }

    fn eat(&mut self, byte: u8) -> Option<()> {
        self.space();
        (self.bytes.get(self.at) == Some(&byte)).then(|| self.at += 1)
    }

    fn value(&mut self) -> Option<Json> {
        self.space();
        match *self.bytes.get(self.at)? {
            b'{' => {
                self.at += 1;
                let mut members = Vec::new();
                if self.eat(b'}').is_some() {
                    return Some(Json::Object(members));
                }
                loop {
                    self.space();
                    let Json::Text(key) = self.value()? else {
                        return None;
                    };
                    self.eat(b':')?;
                    members.push((key, self.value()?));
                    if self.eat(b',').is_none() {
                        self.eat(b'}')?;
                        return Some(Json::Object(members));
                    }
                }
            }
            b'[' => {
                self.at += 1;
                let mut items = Vec::new();
                if self.eat(b']').is_some() {
                    return Some(Json::List(items));
                }
                loop {
                    items.push(self.value()?);
                    if self.eat(b',').is_none() {
                        self.eat(b']')?;
                        return Some(Json::List(items));
                    }
                }
            }
            b'"' => {
                self.at += 1;
                let start = self.at;
                while *self.bytes.get(self.at)? != b'"' {
                    self.at += if self.bytes[self.at] == b'\\' { 2 } else { 1 };
                }
                let text = std::str::from_utf8(&self.bytes[start..self.at]).ok()?;
                self.at += 1;
                Some(Json::Text(text.to_string()))
            }
            b't' => self.word("true", Json::Bool),
            b'f' => self.word("false", Json::Bool),
            b'n' => self.word("null", Json::Null),
            _ => {
                let start = self.at;
                while self
                    .bytes
                    .get(self.at)
                    .is_some_and(|b| b.is_ascii_digit() || b"+-.eE".contains(b))
                {
                    self.at += 1;
                }
                std::str::from_utf8(&self.bytes[start..self.at])
                    .ok()?
                    .parse()
                    .ok()
                    .map(Json::Number)
            }
        }
    }

    fn word(&mut self, word: &str, value: Json) -> Option<Json> {
        self.bytes[self.at..].starts_with(word.as_bytes()).then(|| {
            self.at += word.len();
            value
        })
    }
}
