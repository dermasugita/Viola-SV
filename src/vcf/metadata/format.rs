use anyhow::{Result, anyhow};
use std::collections::HashMap;
use derive_getters::Getters;
use crate::vcf::utils::header_parser::*;
use crate::vcf::metadata::*;

#[derive(Eq, PartialEq, Hash, Clone, Debug)]
pub enum ViolaVcfFormatId {
    PR,
    SR,
    Others(String),
}

impl ViolaVcfFormatId {
    pub fn from_str(s: &str) -> Result<Self> {
        match s {
            "PR" => Ok(ViolaVcfFormatId::PR),
            "SR" => Ok(ViolaVcfFormatId::SR),
            _ => Ok(ViolaVcfFormatId::Others(s.to_string())),
        }
    }
}
#[derive(Getters, Clone)]
pub struct FormatMeta {
    id: ViolaVcfFormatId,
    number: VcfNumber,
    dtype: VcfDataType,
    description: String,
}

impl FormatMeta {
    pub fn new(id: ViolaVcfFormatId, number: VcfNumber, dtype: VcfDataType, description: String) -> Self {
        FormatMeta {
            id,
            number,
            dtype,
            description,
        }
    }
    pub fn from_str(s: &str) -> Result<Self> {
        //##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        // remove ##FORMAT=< and >
        let format_meta = parse_nested_vcf_header(s)?;
        Self::from_hashmap(&format_meta)
    }
    pub fn from_hashmap(h: &HashMap<String, String>) -> Result<Self> {
        let id = match h.get("ID") {
            Some(id) => ViolaVcfFormatId::from_str(id)?,
            None => { return Err(anyhow!("ID is not found in FORMAT. {:?}", h)); }
        };
        let number = match h.get("Number") {
            Some(number) => VcfNumber::from_str(number)?,
            None => { return Err(anyhow!("Number is not found in FORMAT. {:?}", h)); }
        };
        let dtype = match h.get("Type") {
            Some(dtype) => VcfDataType::from_str(dtype)?,
            None => { return Err(anyhow!("Type is not found in FORMAT. {:?}", h)); }
        };
        let description = match h.get("Description") {
            Some(description) => description.to_string(),
            None => { return Err(anyhow!("Description is not found in FORMAT. {:?}", h)); }
        };
        
        Ok(Self::new(
            id,
            number,
            dtype,
            description,
        ))
    }
}