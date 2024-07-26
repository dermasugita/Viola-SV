use anyhow::{Result, anyhow};
use derive_getters::Getters;
use regex::Regex;
use crate::vcf::elements::{DNABaseSeq, Position, End};

#[derive(Debug, Clone, PartialEq, Eq, Getters)]
pub struct Breakend {
    t: DNABaseSeq,
    this_end: End,
    other_end_pos: Option<Position>,
    other_end_extend_to: Option<End>, // None means DEADEND
}

impl Breakend {
    pub fn from_str(s: &str) -> Result<Self> {
        let re = Regex::new(r"[\[\]]")?;
        // check s contains at least one of [ or ]
        if !re.is_match(s) {
            if s.starts_with(".") {
                return Ok(Self {
                    t: DNABaseSeq::from_str(&s[1..])?,
                    this_end: End::FivePrime,
                    other_end_pos: None,
                    other_end_extend_to: None
                })
            }
            if s.contains(".") {
                return Ok(Self {
                    t: DNABaseSeq::from_str(&s[..s.len()-1])?,
                    this_end: End::ThreePrime,
                    other_end_pos: None,
                    other_end_extend_to: None
                })
            }
            return Err(anyhow!("Invalid Breakend format."))
        }
        if s.starts_with("]") {
            let mut it = s[1..].split(']');
            let other_end_pos = Position::from_colon_separated(it.next().unwrap());
            let t = DNABaseSeq::from_str(it.next().unwrap())?;
            return Ok(Self {
                t,
                this_end: End::FivePrime,
                other_end_pos: Some(other_end_pos),
                other_end_extend_to: Some(End::FivePrime)
            })
        }
        if s.starts_with("[") {
            let mut it = s[1..].split('[');
            let other_end_pos = Position::from_colon_separated(it.next().unwrap());
            let t = DNABaseSeq::from_str(it.next().unwrap())?;
            return Ok(Self {
                t,
                this_end: End::FivePrime,
                other_end_pos: Some(other_end_pos),
                other_end_extend_to: Some(End::ThreePrime)
            })
        }
        if s.contains("]") {
            let mut it = s.split(']');
            let t = DNABaseSeq::from_str(it.next().unwrap())?;
            let other_end_pos = Position::from_colon_separated(it.next().unwrap());
            return Ok(Self {
                t,
                this_end: End::ThreePrime,
                other_end_pos: Some(other_end_pos),
                other_end_extend_to: Some(End::FivePrime)
            })
        }
        let mut it = s.split('[');
        let t = DNABaseSeq::from_str(it.next().unwrap())?;
        let other_end_pos = Position::from_colon_separated(it.next().unwrap());
        Ok(Self {
            t,
            this_end: End::ThreePrime,
            other_end_pos: Some(other_end_pos),
            other_end_extend_to: Some(End::ThreePrime)
        })
    }
    pub fn to_string(&self) -> String {
        if let None = self.other_end_pos {
            match self.this_end() {
                End::FivePrime => return format!(".{}", self.t().to_string()),
                End::ThreePrime => return format!("{}.", self.t().to_string()),
            }
        }
        let other_end_pos = self.other_end_pos().as_ref().unwrap();
        let other_end_extend_to = self.other_end_extend_to().as_ref().unwrap();
        match self.this_end() {
            End::FivePrime => {
                match other_end_extend_to {
                    End::FivePrime => format!("]{}]{}", other_end_pos.get_colon_separated(), self.t().to_string()),
                    End::ThreePrime => format!("[{}[{}", other_end_pos.get_colon_separated(), self.t().to_string()),
                }
            },
            End::ThreePrime => {
                match other_end_extend_to {
                    End::FivePrime => format!("{}]{}[", self.t().to_string(), other_end_pos.get_colon_separated()),
                    End::ThreePrime => format!("{}[{}[", self.t().to_string(), other_end_pos.get_colon_separated()),
                }
            }

        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AltElement {
    BaseSeq(DNABaseSeq),
    Asterisk,
    Missing,
    AngleBracket(String),
    Breakend(Breakend),
}

impl AltElement {
    pub fn from_str(s: &str) -> Result<Self> {
        if s == "*" {
            return Ok(Self::Asterisk)
        }
        if s == "." || s == "MISSING" {
            return Ok(Self::Missing)
        }
        if s.starts_with("<") && s.ends_with(">") {
            return Ok(Self::AngleBracket(s[1..s.len()-1].to_string()))
        }
        if s.contains("[") || s.contains("]") || s.contains(".") {
            return Ok(Self::Breakend(Breakend::from_str(s)?))
        }
        Ok(Self::BaseSeq(DNABaseSeq::from_str(s)?))
    }
    pub fn to_string(&self) -> String {
        match self {
            Self::BaseSeq(base_seq) => base_seq.to_string(),
            Self::Asterisk => "*".to_string(),
            Self::Missing => ".".to_string(),
            Self::AngleBracket(s) => format!("<{}>", s),
            Self::Breakend(breakend) => breakend.to_string(),
        }
    }
}