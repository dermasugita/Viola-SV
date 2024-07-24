use std::collections::HashMap;
use crate::vcf::metadata::*;
use crate::vcf::utils::header_parser::parse_nested_vcf_header;
use anyhow::{Result, anyhow};
use url::Url;

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

impl VcfNumber {
    pub fn from_str(s: &str) -> Result<Self> {
        match s {
            "A" => Ok(VcfNumber::A),
            "G" => Ok(VcfNumber::G),
            "R" => Ok(VcfNumber::R),
            "." => Ok(VcfNumber::Dot),
            _ => {
                match s.parse::<i32>() {
                    Ok(n) => Ok(VcfNumber::Number(n)),
                    Err(_) => Err(anyhow!("Unknown VcfNumber has been found. {}", s)),
                }
            }
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

impl InfoMeta {
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

impl FormatMeta {
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

impl AltMeta {
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