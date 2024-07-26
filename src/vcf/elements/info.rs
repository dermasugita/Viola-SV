use std::collections::HashMap;
use anyhow::{Result, anyhow};
use crate::vcf::elements::data::VcfData;
use crate::vcf::metadata::*;

#[derive(Debug, Clone, PartialEq, Eq)]
struct TypeNaiveInfo (HashMap<ViolaVcfInfoId, Option<Vec<String>>>);

#[derive(Debug, Clone)]
pub struct Info (HashMap<ViolaVcfInfoId, Vec<VcfData>>);
impl TypeNaiveInfo {
    // "IMPRECISE;SVTYPE=DEL;SVLEN=-100000;END=200000;CIPOS=-51,52;CIEND=-51,52;SOMATIC;SOMATICSCORE=10"
    pub fn from_str(s: &str) -> Result<Self> {
        let mut info: HashMap<ViolaVcfInfoId, Option<Vec<String>>> = HashMap::new();
        for info_str in s.split(";") {
            let mut kv_it = info_str.split("=");
            let key = match kv_it.next() {
                Some(k) => {
                    ViolaVcfInfoId::from_str(k)?
                },
                None => { return Err(anyhow!("Empty INFO has been found. {}", info_str)); }
            };
            let value = match kv_it.next() {
                Some(v) => {
                    Some(v.split(',').map(|s| s.to_string()).collect::<Vec<String>>())
                },
                None => None,
            };
            info.insert(key, value); 
        }
        Ok(Self(info))
    }
    //fn to_type_aware_info(&self, vcf_info_meta: &[InfoMeta]) -> Result<Info> {
    pub fn to_type_aware_info(&self, vcf_info_meta: &[InfoMeta]) -> Result<Info> {
        let mut info_meta_map: HashMap<ViolaVcfInfoId, InfoMeta> = HashMap::new();
        for meta in vcf_info_meta {
            info_meta_map.insert(meta.id().clone(), meta.clone());
        }
        let mut info: HashMap<ViolaVcfInfoId, Vec<VcfData>> = HashMap::new();
        for (k, v) in &self.0 {
            let meta = match info_meta_map.get(k) {
                Some(meta) => meta,
                None => { return Err(anyhow!("No meta data found for INFO: {:?}", k)); }
            };
            let mut data: Vec<VcfData> = Vec::new();
            match v {
                Some(v) => {
                    for s in v {
                        data.push(VcfData::from_str(&Some(s.to_string()), &meta.dtype())?);
                    }
                },
                None => {
                    data.push(VcfData::from_str(&None, &meta.dtype())?);
                }
            }
            info.insert(k.clone(), data);
        }

        Ok(Info(info))
    }
}

impl Info {
    pub fn from_str(s: &str, vcf_info_meta: &[InfoMeta]) -> Result<Self> {
        let type_naive_info = TypeNaiveInfo::from_str(s)?;
        type_naive_info.to_type_aware_info(vcf_info_meta)
    }
    pub fn get(&self, key: &str) -> Option<&Vec<VcfData>> {
        self.0.get(&ViolaVcfInfoId::from_str(key).unwrap())
    }
    pub fn get_hashmap(&self) -> &HashMap<ViolaVcfInfoId, Vec<VcfData>> {
        &self.0
    }
}