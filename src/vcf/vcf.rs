use chrono::{NaiveDate};
use anyhow::{Result, anyhow};
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::path::PathBuf;
use pyo3::prelude::*;
use crate::vcf::metadata::*;
use crate::vcf::records::*;
use crate::vcf::utils::header_parser::{parse_vcf_header, ABCVcfHeader};

//#[derive(FromPyObject)]
pub struct Vcf<T> where T: VcfRecords {
    fileformat: String,
    filedate: Option<NaiveDate>,
    source: Option<String>,
    reference: Option<String>,
    phasing: Option<String>,
    contig: Option<Vec<ContigMeta>>,
    kv: Option<Vec<VcfKV>>,
    nested_kv: Option<Vec<VcfNestedKV>>,
    info: Vec<InfoMeta>,
    filter: Vec<FilterMeta>,
    format: Vec<FormatMeta>,
    alt: Vec<AltMeta>,
    meta: Vec<Meta>,
    sample: Vec<SampleMeta>,
    required_header: Vec<VcfColumnRequiredHeader>,
    optional_header: Option<Vec<VcfColumnOptionalHeader>>,
    records: T,
}

impl<T> Vcf<T> where T: VcfRecords {
    pub fn read_from_file(file: &PathBuf) -> Result<Self> {
        let mut fileformat: String = "".to_string();
        let mut filedate: Option<NaiveDate> = None;
        let mut source: Option<String> = None;
        let mut reference: Option<String> = None;
        let mut phasing: Option<String> = None;
        let mut contig: Option<Vec<ContigMeta>> = None;
        let mut kv: Option<Vec<VcfKV>> = None;
        let mut nested_kv: Option<Vec<VcfNestedKV>> = None;
        let mut info: Vec<InfoMeta> = Vec::new();
        let mut filter: Vec<FilterMeta> = Vec::new();
        let mut format: Vec<FormatMeta> = Vec::new();
        let mut alt: Vec<AltMeta> = Vec::new();
        let mut meta: Vec<Meta> = Vec::new();
        let mut sample: Vec<SampleMeta> = Vec::new();
        let mut required_header: Vec<VcfColumnRequiredHeader> = Vec::new();
        let mut optional_header: Option<Vec<VcfColumnOptionalHeader>> = None;

        for line in BufReader::new(File::open(file)?).lines() {
            let line = line?;
            if !(&line).starts_with("#") { break; }
            let parsed_line = parse_vcf_header(&line)?;
            match parsed_line {
                ABCVcfHeader::NestedKV(nkv) => {
                    match nkv.key().as_str() {
                        "INFO" => {
                            let info_meta = InfoMeta::from_hashmap(nkv.value())?;
                            info.push(info_meta);
                            // drop INFO from nkv
                        },
                        "FILTER" => {
                            let filter_meta = FilterMeta::from_hashmap(nkv.value())?;
                            filter.push(filter_meta);
                        },
                        "FORMAT" => {
                            let format_meta = FormatMeta::from_hashmap(nkv.value())?;
                            format.push(format_meta);
                        },
                        "ALT" => {
                            let alt_meta = AltMeta::from_hashmap(nkv.value())?;
                            alt.push(alt_meta);
                        },
                        "META" => {
                            let meta_meta = Meta::from_hashmap(nkv.value())?;
                            meta.push(meta_meta);
                        },
                        "CONTIG" => {
                            let contig_meta = ContigMeta::from_hashmap(nkv.value())?;
                            match &mut contig {
                                Some(some_contig) => {
                                    some_contig.push(contig_meta);
                                },
                                None => {
                                    contig = Some(vec![contig_meta]);
                                }
                            }
                        },
                        "SAMPLE" => {
                            let mut meta_ids = Vec::new();
                            if meta.len() == 0 {
                                return Err(anyhow!("META should be defined before SAMPLE"));
                            }
                            for m in &meta {
                                meta_ids.push(m.id().clone());
                            }
                            let sample_meta = SampleMeta::from_hashmap(nkv.value(), &meta_ids)?;
                            sample.push(sample_meta);
                        },
                        _ => {
                            match &mut nested_kv {
                                Some(some_nested_kv) => {
                                    some_nested_kv.push(nkv);
                                },
                                None => {
                                    nested_kv = Some(vec![nkv]);
                                }
                            }
                        }
                    }
                },
                ABCVcfHeader::KV(kv_) => {
                    match kv_.key().as_str() {
                        "fileformat" => {
                            fileformat = kv_.value().to_string();
                        },
                        "filedate" => {
                            filedate = Some(NaiveDate::parse_from_str(kv_.value(), "%Y%m%d")?);
                        },
                        "source" => {
                            source = Some(kv_.value().to_string());
                        },
                        "reference" => {
                            reference = Some(kv_.value().to_string());
                        },
                        "phasing" => {
                            phasing = Some(kv_.value().to_string());
                        },
                        _ => {
                            match &mut kv {
                                Some(some_kv) => {
                                    some_kv.push(kv_);
                                },
                                None => {
                                    kv = Some(vec![kv_]);
                                }
                            }
                        }
                    }
                }
            }
        }
       let records = T::from_path(&PathBuf::from(file))?; 
       Ok(
        Self {
            fileformat,
            filedate,
            source,
            reference,
            phasing,
            contig,
            kv,
            nested_kv,
            info,
            filter,
            format,
            alt,
            meta,
            sample,
            required_header,
            optional_header,
            records,
        }
       )
    }
    pub fn get_by_id(&self, id: &str) -> Option<VcfRecord> {
        self.records.get_by_id(id)
    }
}