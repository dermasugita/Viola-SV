use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::path::PathBuf;
use std::fs::File;
mod vcf;
mod vcf_impls;
use vcf::vcf::Vcf;
use vcf::mrec::InMemoryVcfRecords;
use vcf_impls::VcfInMemory;

#[pyfunction]
fn rust_bounded(invalue: &str) -> PyResult<String> {
    Ok(invalue.to_string())
}

/// Formats the sum of two numbers as string.
#[pyfunction]
fn read_vcf(path: PathBuf) -> PyResult<VcfInMemory> {
    VcfInMemory::read_from_file(path)
}

/// A Python module implemented in Rust.
#[pymodule]
fn rs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rust_bounded, m)?)?;
    m.add_function(wrap_pyfunction!(read_vcf, m)?)?;
    Ok(())
}
