use std::path::Path;
use anyhow::Result;
use derive_getters::Getters;
use pyo3::prelude::*;
use crate::vcf::elements::*;

#[pyclass]
#[derive(Getters, Clone)]
pub struct VcfRecord {
    chrom: Contig,
    pos: i128,
    id: String,
    ref_: DNABaseSeq,
    alt: Vec<AltElement>,
    qual: Option<f32>,
    filter: Vec<Filter>,
    info: Info,
    format: Format,
}

impl VcfRecord {
    pub fn new(chrom: Contig, pos: i128, id: String, ref_: DNABaseSeq, alt: Vec<AltElement>, qual: Option<f32>, filter: Vec<Filter>, info: Info, format: Format) -> Self {
        Self {
            chrom,
            pos,
            id,
            ref_,
            alt,
            qual,
            filter,
            info,
            format,
        }
    }
}
#[pymethods]
impl VcfRecord {
    pub fn __repr__(&self) -> PyResult<String> {
        Ok(format!("VcfRecord(chrom={:?}, pos={}, id={}, ref_={:?}, alt={:?}, qual={:?}, filter={:?}, info={:?}, format={:?})", self.chrom, self.pos, self.id, self.ref_, self.alt, self.qual, self.filter, self.info, self.format))
    }
}


pub trait VcfRecords: Sized {
    fn from_path<P: AsRef<Path>>(path: &P) -> Result<Self>;
    fn get_by_id(&self, id: &str) -> Option<VcfRecord>;
}
