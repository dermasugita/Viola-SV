
use std::collections::HashMap;
use anyhow::{Result, anyhow};
use derive_getters::Getters;
use crate::vcf::utils::header_parser::*;
use crate::vcf::metadata::*;
// Below 3 structs are used for sample meta data
#[derive(Clone, Getters)]
pub struct MetaId {
    name: String,
}
impl MetaId {
    pub fn new(name: String) -> Self {
        MetaId {
            name,
        }
    }
}
#[derive(Clone, Getters)]
pub struct Meta {
    id: MetaId,
    dtype: VcfDataType,
    number: VcfNumber,
    values: Vec<String>,
    description: String,
}
impl Meta {
    pub fn new(id: MetaId, dtype: VcfDataType, number: VcfNumber, values: Vec<String>, description: String) -> Self {
        Meta {
            id,
            dtype,
            number,
            values,
            description,
        }
    }
}

pub struct SampleMeta {
    id: String,
    description: String,
    kv: Vec<VcfKV<MetaId, String>>,
}
impl SampleMeta {
    pub fn new(id: String, description: String, kv: Vec<VcfKV<MetaId, String>>) -> Self {
        SampleMeta {
            id,
            description,
            kv,
        }
    }
}

impl Meta {
    pub fn from_str(s: &str) -> Result<Self> {
        let meta = parse_nested_vcf_header(s)?;
        Self::from_hashmap(&meta)
    }
    pub fn from_hashmap(h: &HashMap<String, String>) -> Result<Self> {
        let id = match h.get("ID") {
            Some(id) => MetaId::new(id.to_string()),
            None => { return Err(anyhow!("ID is not found in META. {:?}", h)); }
        };
        let description = match h.get("Description") {
            Some(description) => description.to_string(),
            None => { return Err(anyhow!("Description is not found in META. {:?}", h)); }
        };
        let number = match h.get("Number") {
            Some(number) => VcfNumber::from_str(number)?,
            None => { return Err(anyhow!("Number is not found in META. {:?}", h)); }
        };
        let dtype = match h.get("Type") {
            Some(dtype) => VcfDataType::from_str(dtype)?,
            None => { return Err(anyhow!("Type is not found in META. {:?}", h)); }
        };
        let values = match h.get("Values") {
            Some(values) => {
                let values = values.trim_start_matches("[").trim_end_matches("]").replace(" ", "");
                values.split(",").map(|s| s.to_string()).collect::<Vec<String>>()},
            None => { return Err(anyhow!("Values is not found in META. {:?}", h)); }
        };
        
        Ok(Self::new(
            id,
            dtype,
            number,
            values,
            description,
        ))
    }
}

impl SampleMeta {
    pub fn from_str(s: &str, meta_ids: &[MetaId]) -> Result<Self> {
        let sample_meta = parse_nested_vcf_header(s)?;
        Self::from_hashmap(&sample_meta, meta_ids)
    }
    pub fn from_hashmap(s: &HashMap<String, String>, meta_ids: &[MetaId]) -> Result<Self> {
        let id = match s.get("ID") {
            Some(id) => id.to_string(),
            None => { return Err(anyhow!("ID is not found in SAMPLE. {:?}", s)); }
        };
        let description = match s.get("Description") {
            Some(description) => description.to_string(),
            None => { return Err(anyhow!("Description is not found in SAMPLE. {:?}", s)); }
        };
        let mut kv: Vec<VcfKV<MetaId, String>> = Vec::new();
        for meta_id in meta_ids {
            let val = match s.get(meta_id.name()) {
                Some(kv) => kv.to_string(),
                None => { return Err(anyhow!("Meta data is not found in SAMPLE. {:?}", s)); }
            };
            kv.push(VcfKV::new(meta_id.clone(), val));
        }
        
        Ok(Self::new(
            id,
            description,
            kv,
        ))
    }
}