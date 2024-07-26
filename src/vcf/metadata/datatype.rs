use anyhow::{Result, anyhow};
#[derive(Copy, Clone)]
pub enum VcfDataType {
    Integer,
    Float,
    Flag,
    Character,
    String,
}

impl VcfDataType {
    pub fn from_str(s: &str) -> Result<Self> {
        match s {
            "Flag" => Ok(VcfDataType::Flag),
            "Integer" => Ok(VcfDataType::Integer),
            "Float" => Ok(VcfDataType::Float),
            "Character" => Ok(VcfDataType::Character),
            "String" => Ok(VcfDataType::String),
            _ => Err(anyhow!("Unknown VcfDataType has been found. {}", s)),
        }
    }
}