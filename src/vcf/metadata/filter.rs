use std::collections::HashMap;
use anyhow::{Result, anyhow};
use derive_getters::Getters;
use crate::vcf::utils::header_parser::*;


#[derive(Getters, Clone)]
pub struct FilterMeta {
    id: String,
    description: String,
}
impl FilterMeta {
    pub fn new(id: String, description: String) -> Self {
        FilterMeta {
            id,
            description,
        }
    }
}
impl FilterMeta {
    pub fn from_str(s: &str) -> Result<Self> {
        //##FILTER=<ID=q10,Description="Quality below 10">
        // remove ##FILTER=< and >
        let filter_meta = parse_nested_vcf_header(s)?;
        Self::from_hashmap(&filter_meta)
    }
    pub fn from_hashmap(h: &HashMap<String, String>) -> Result<Self> {
        let id = match h.get("ID") {
            Some(id) => id.to_string(),
            None => { return Err(anyhow!("ID is not found in FILTER. {:?}", h)); }
        };
        let description = match h.get("Description") {
            Some(description) => description.to_string(),
            None => { return Err(anyhow!("Description is not found in FILTER. {:?}", h)); }
        };
        
        Ok(Self::new(
            id,
            description,
        ))
    }
}