use std::path::PathBuf;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use crate::vcf::Vcf;
use crate::vcf::records::*;
use crate::vcf::callers::*;

#[pyclass]
pub struct VcfInMemory {
    inner: Vcf<InMemoryVcfRecords>,
}

#[pymethods]
impl VcfInMemory {
    #[new]
    pub fn read_from_file(path: PathBuf, caller: &str) -> PyResult<Self> {
        match Vcf::<InMemoryVcfRecords>::read_from_file(&path, caller) {
            Ok(vcf) => Ok(Self { inner: vcf }),
            Err(e) => Err(PyValueError::new_err(e.to_string())),
        }
    }
    pub fn sv_count(&self) -> usize {
        self.inner.sv_count()
    }
    pub fn contigs(&self) -> Vec<String> {
        self.inner.contigs()
    }
    pub fn ids(&mut self) -> Vec<String> {
        self.inner.ids()
    }
    pub fn get_positions_table(&mut self) -> PositionsTable {
        self.inner.get_positions_table().clone()
    }
    pub fn get_info_tables(&mut self) -> InfoTables {
        self.inner.get_info_tables().clone()
    }
    pub fn get_format_tables(&mut self) -> FormatTables {
        self.inner.get_format_tables().clone()
    }
}