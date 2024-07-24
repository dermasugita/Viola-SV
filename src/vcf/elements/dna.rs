use anyhow::{Result, anyhow};
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum IUPACAmbiguityCodes {
    R,
    Y,
    S,
    W,
    K,
    M,
    B,
    D,
    H,
    V,
    GAP,
}
impl IUPACAmbiguityCodes {
    pub fn from_char(c: char) -> Result<Self> {
        let c = c.to_uppercase().next().unwrap();
        match c {
            'R' => Ok(Self::R),
            'Y' => Ok(Self::Y),
            'S' => Ok(Self::S),
            'W' => Ok(Self::W),
            'K' => Ok(Self::K),
            'M' => Ok(Self::M),
            'B' => Ok(Self::B),
            'D' => Ok(Self::D),
            'H' => Ok(Self::H),
            'V' => Ok(Self::V),
            '-' => Ok(Self::GAP),
            '.' => Ok(Self::GAP),
            _ => Err(anyhow!("Invalid IUPAC ambiguity code character."))
        }
    }
    pub fn to_char(&self) -> char {
        match self {
            Self::R => 'R',
            Self::Y => 'Y',
            Self::S => 'S',
            Self::W => 'W',
            Self::K => 'K',
            Self::M => 'M',
            Self::B => 'B',
            Self::D => 'D',
            Self::H => 'H',
            Self::V => 'V',
            Self::GAP => '-',
        }
    }
}
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum DNABase {
    A,
    T,
    G,
    C,
    N,
    IUPAC(IUPACAmbiguityCodes),
}

impl DNABase {
    pub fn from_char(c: char) -> Result<Self> {
        let c = c.to_uppercase().next().unwrap();
        match c {
            'A' => Ok(Self::A),
            'T' => Ok(Self::T),
            'G' => Ok(Self::G),
            'C' => Ok(Self::C),
            'N' => Ok(Self::N),
            _ => {
                match IUPACAmbiguityCodes::from_char(c) {
                    Ok(iupac) => Ok(Self::IUPAC(iupac)),
                    Err(e) => Err(e),
                }
            }
        }
    }
    pub fn to_char(&self) -> char {
        match self {
            Self::A => 'A',
            Self::T => 'T',
            Self::G => 'G',
            Self::C => 'C',
            Self::N => 'N',
            Self::IUPAC(iupac) => iupac.to_char(),
        }
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct DNABaseSeq(Vec<DNABase>);

impl DNABaseSeq {
    pub fn from_str(s: &str) -> Result<Self> {
        let mut seq = Vec::new();
        for c in s.chars() {
            seq.push(DNABase::from_char(c)?);
        }
        Ok(Self(seq))
    }
    pub fn to_string(&self) -> String {
        self.0.iter().map(|b| b.to_char()).collect::<String>()
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub enum End {
    ThreePrime,
    FivePrime,
}