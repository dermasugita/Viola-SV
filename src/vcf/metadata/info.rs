use anyhow::{Result, anyhow};
use std::collections::HashMap;
use crate::vcf::metadata::*;
use crate::vcf::utils::header_parser::*;
use derive_getters::Getters;

#[derive(Eq, PartialEq, Hash, Clone, Debug)]
pub enum ViolaVcfInfoId {
    CIPOS,
    CIEND,
    IMPRECISE,
    PRECISE,
    SVTYPE,
    Others(String),
}
impl ViolaVcfInfoId {
    pub fn from_str(s: &str) -> Result<Self> {
        match s {
            "CIPOS" => Ok(ViolaVcfInfoId::CIPOS),
            "CIEND" => Ok(ViolaVcfInfoId::CIEND),
            "IMPRECISE" => Ok(ViolaVcfInfoId::IMPRECISE),
            "PRECISE" => Ok(ViolaVcfInfoId::PRECISE),
            "SVTYPE" => Ok(ViolaVcfInfoId::SVTYPE),
            _ => Ok(ViolaVcfInfoId::Others(s.to_string())),
        }
    }
    pub fn to_string(&self) -> String {
        match self {
            ViolaVcfInfoId::CIPOS => "CIPOS".to_string(),
            ViolaVcfInfoId::CIEND => "CIEND".to_string(),
            ViolaVcfInfoId::IMPRECISE => "IMPRECISE".to_string(),
            ViolaVcfInfoId::PRECISE => "PRECISE".to_string(),
            ViolaVcfInfoId::SVTYPE => "SVTYPE".to_string(),
            ViolaVcfInfoId::Others(s) => s.to_string(),
        }
    }
}

#[derive(Getters, Clone)]
pub struct InfoMeta {
    id: ViolaVcfInfoId,
    number: VcfNumber,
    dtype: VcfDataType,
    description: String,
    source: Option<String>,
    version: Option<String>,
}
impl InfoMeta {
    pub fn new(id: ViolaVcfInfoId, number: VcfNumber, dtype: VcfDataType, description: String, source: Option<String>, version: Option<String>) -> Self {
        InfoMeta {
            id,
            number,
            dtype,
            description,
            source,
            version,
        }
    }
    pub fn from_str(s: &str) -> Result<Self> {
        //##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
        // remove ##INFO=< and >
        let info_meta = parse_nested_vcf_header(s)?;
        Self::from_hashmap(&info_meta)
    }
    pub fn from_hashmap(h: &HashMap<String, String>) -> Result<Self> {
        let id = match h.get("ID") {
            Some(id) => ViolaVcfInfoId::from_str(id)?,
            None => { return Err(anyhow!("ID is not found in INFO. {:?}", h)); }
        };
        let number = match h.get("Number") {
            Some(number) => VcfNumber::from_str(number)?,
            None => { return Err(anyhow!("Number is not found in INFO. {:?}", h)); }
        };
        let dtype = match h.get("Type") {
            Some(dtype) => VcfDataType::from_str(dtype)?,
            None => { return Err(anyhow!("Type is not found in INFO. {:?}", h)); }
        };
        let description = match h.get("Description") {
            Some(description) => description.to_string(),
            None => { return Err(anyhow!("Description is not found in INFO. {:?}", h)); }
        };
        let source = h.get("Source").map(|s| s.to_string());
        let version = h.get("Version").map(|s| s.to_string());
        
        Ok(Self::new(
            id,
            number,
            dtype,
            description,
            source,
            version,
        ))
    }

}