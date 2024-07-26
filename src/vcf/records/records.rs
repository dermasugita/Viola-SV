use std::path::Path;
use anyhow::Result;
use derive_getters::Getters;
use pyo3::prelude::*;
use crate::vcf::{callers::*, elements::*};

#[pyclass]
#[derive(Getters, Clone)]
pub struct VcfRecord {
    chrom: Contig,
    pos: i64,
    id: String,
    ref_: DNABaseSeq,
    alt: Vec<AltElement>,
    qual: Option<f32>,
    filter: Vec<Filter>,
    info: Info,
    format: Format,
}

impl VcfRecord {
    pub fn new(chrom: Contig, pos: i64, id: String, ref_: DNABaseSeq, alt: Vec<AltElement>, qual: Option<f32>, filter: Vec<Filter>, info: Info, format: Format) -> Self {
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


pub trait VcfRecords: Sized + Iterator<Item=VcfRecord> {
    fn from_path<P: AsRef<Path>>(path: &P, caller: &SupportedCallers) -> Result<Self>;
    fn get_by_ids(&self, ids: &[&str]) -> Result<Vec<VcfRecord>>;
    fn remove_by_id(&mut self, id: &str) -> Option<VcfRecord>;
    fn count(&self) -> usize;
    fn get_positions_table(&self) -> &PositionsTable;
    fn get_info_tables(&self) -> &InfoTables;
    fn get_format_tables(&self) -> &FormatTables;
}