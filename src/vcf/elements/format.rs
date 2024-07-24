use std::{collections::HashMap, path::Path};
use std::path::PathBuf;
use anyhow::{Result, anyhow};
use derive_getters::Getters;
use regex::Regex;
use crate::vcf::metadata::*;
use crate::vcf::elements::data::VcfData;
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TypeNaiveFormatUnit(HashMap<ViolaVcfFormatId, Vec<String>>);
#[derive(Debug, Clone)]
pub struct FormatUnit(HashMap<ViolaVcfFormatId, Vec<VcfData>>);
#[derive(Debug, Clone)]
pub struct Format(HashMap<String, FormatUnit>); // sample id -> format id -> format values

impl TypeNaiveFormatUnit {
    pub fn from_str(guide: &str, values: &str) -> Result<Self> { // guide: "PR:SR", values: "1:2"
        let mut format_unit: HashMap<ViolaVcfFormatId, Vec<String>> = HashMap::new();
        for (guide, value) in guide.split(":").zip(values.split(":")) {
            let key = ViolaVcfFormatId::from_str(guide)?;
            let value = value.split(",").map(|s| s.to_string()).collect::<Vec<String>>();
            format_unit.insert(key, value);
        }
        Ok(Self(format_unit))
    }
    pub fn to_type_aware_format_unit(&self, vcf_format_meta: &[FormatMeta]) -> Result<FormatUnit> {
        let mut format_meta_map: HashMap<ViolaVcfFormatId, FormatMeta> = HashMap::new();
        for meta in vcf_format_meta {
            format_meta_map.insert(meta.id().clone(), meta.clone());
        }
        let mut format_unit: HashMap<ViolaVcfFormatId, Vec<VcfData>> = HashMap::new();
        for (k, v) in &self.0 {
            let meta = match format_meta_map.get(k) {
                Some(meta) => meta,
                None => { return Err(anyhow!("No meta data found for FORMAT: {:?}", k)); }
            };
            let mut data: Vec<VcfData> = Vec::new();
            for s in v {
                data.push(VcfData::from_str(&Some(s.to_string()), &meta.dtype())?);
            }
            format_unit.insert(k.clone(), data);
        }
        Ok(FormatUnit(format_unit))
    }
}

impl FormatUnit {
    pub fn from_str(guids: &str, values: &str, vcf_format_meta: &[FormatMeta]) -> Result<Self> {
        let type_naive_format_unit = TypeNaiveFormatUnit::from_str(guids, values)?;
        type_naive_format_unit.to_type_aware_format_unit(vcf_format_meta)
    }
    pub fn get(&self, format_id: &str) -> Option<&Vec<VcfData>> {
        self.0.get(&ViolaVcfFormatId::from_str(format_id).unwrap())
    }

}

impl Format {
    pub fn new() -> Self {
        Self(HashMap::new())
    }
    pub fn insert(&mut self, sample_id: String, format_unit: FormatUnit) {
        self.0.insert(sample_id, format_unit);
    }
    pub fn get(&self, sample_id: &str) -> Option<&FormatUnit> {
        self.0.get(sample_id)
    }
}