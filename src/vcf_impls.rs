use chrono::{NaiveDate};
use pyo3::exceptions::PyValueError;
use anyhow::{Result, anyhow};
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::path::PathBuf;
use pyo3::prelude::*;
use crate::vcf::vcf::Vcf;
use crate::vcf::metadata::*;
use crate::vcf::records::*;
use crate::vcf::utils::header_parser::{parse_vcf_header, ABCVcfHeader};
use crate::vcf::mrec::InMemoryVcfRecords;

macro_rules! create_vcf {
    ($name: ident, $type: ident) => {
        #[pyclass]
        pub struct $name {
            inner: Vcf<$type>,
        }
        #[pymethods]
        impl $name {
            #[new]
            pub fn read_from_file(path: PathBuf) -> PyResult<Self> {
                match Vcf::<$type>::read_from_file(&path) {
                    Ok(vcf) => Ok(Self { inner: vcf }),
                    Err(e) => Err(PyValueError::new_err(e.to_string())),
                }
            }
            pub fn get_by_id(&self, id: &str) -> Option<VcfRecord> {
                self.inner.get_by_id(id)
            }
        }
    };
}


create_vcf!(VcfInMemory, InMemoryVcfRecords);