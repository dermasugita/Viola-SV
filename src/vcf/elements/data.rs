
use anyhow::{Result, anyhow};
use regex::Regex;
use crate::vcf::metadata::*;
#[derive(Clone, Debug, PartialEq)]
pub enum VcfData {
    Integer(i32),
    Float(Option<f32>),
    Flag(bool),
    Character(char),
    String(String),
}
impl VcfData {
    pub fn from_str(s: &Option<String>, dtype: &VcfDataType) -> Result<Self> {
        match s {
            Some(s) => {
                match dtype {
                    VcfDataType::Integer => {
                        Ok(Self::Integer(s.parse::<i32>()?))
                    },
                    VcfDataType::Float => {
                        let positive_inf_re = Regex::new(r"^+?(INF|INFINITY)")?;
                        let negative_inf_re = Regex::new(r"^-(INF|INFINITY)")?;
                        if positive_inf_re.is_match(s) {
                            return Ok(Self::Float(Some(f32::INFINITY)))
                        }
                        if negative_inf_re.is_match(s) {
                            return Ok(Self::Float(Some(f32::NEG_INFINITY)))
                        }
                        if s == "NAN" {
                            return Ok(Self::Float(None))
                        }
                        Ok(Self::Float(Some(s.parse::<f32>()?)))
                    },
                    VcfDataType::Flag => {
                        Err(anyhow!("Flag type should not have value."))
                    },
                    VcfDataType::Character => {
                        Ok(Self::Character(s.chars().next().unwrap()))
                    },
                    VcfDataType::String => {
                        Ok(Self::String(s.to_string()))
                    },
                }
            },
            None => {
                if let VcfDataType::Flag = dtype {
                    Ok(Self::Flag(true))
                } else {
                    Err(anyhow!("Empty value has been found."))
                }
            }
        }
    }
}