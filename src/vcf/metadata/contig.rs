use anyhow::{Result, anyhow};
use url::Url;
use std::collections::HashMap;
use derive_getters::Getters;
use crate::vcf::utils::header_parser::*;

#[derive(Getters, Clone)]
pub struct ContigMeta {
    id: String,
    length: Option<i128>,
    assembly: Option<String>,
    url: Option<Url>,
    md5: Option<String>,
    species: Option<String>,
    taxonomy: Option<String>,
}
impl ContigMeta {
    pub fn new(id: String, length: Option<i128>, assembly: Option<String>, url: Option<Url>, md5: Option<String>, species: Option<String>, taxonomy: Option<String>) -> Self {
        ContigMeta {
            id,
            length,
            assembly,
            url,
            md5,
            species,
            taxonomy,
        }
    }

}
impl ContigMeta {
    pub fn from_str(s: &str) -> Result<Self> {
        let contig_meta = parse_nested_vcf_header(s)?;
        Self::from_hashmap(&contig_meta)
    }
    pub fn from_hashmap(h: &HashMap<String, String>) -> Result<Self> {
        let id = match h.get("ID") {
            Some(id) => id.to_string(),
            None => { return Err(anyhow!("ID is not found in contig. {:?}", h)); }
        };
        let length = h.get("length").map(|s| s.parse::<i128>().unwrap());
        let assembly = h.get("assembly").map(|s| s.to_string());
        let url = h.get("url").map(|s| Url::parse(s).unwrap());
        let md5 = h.get("md5").map(|s| s.to_string());
        let species = h.get("species").map(|s| s.to_string());
        let taxonomy = h.get("taxonomy").map(|s| s.to_string());
        
        Ok(Self::new(
            id,
            length,
            assembly,
            url,
            md5,
            species,
            taxonomy,
        ))
    }
}