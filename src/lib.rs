use pyo3::prelude::*;
use std::fs::File;
mod vcf;
use vcf::vcf::Vcf;
use vcf::mrec::InMemoryVcfRecords;

#[pyfunction]
fn rust_bounded(invalue: &str) -> PyResult<String> {
    Ok(invalue.to_string())
}

/// Formats the sum of two numbers as string.
#[pyfunction]
fn read_vcf(path: &str) -> PyResult<()> {
    Vcf::<InMemoryVcfRecords>::read_from_file(path)?;
    Ok(())
}

/// A Python module implemented in Rust.
#[pymodule]
fn rs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rust_bounded, m)?)?;
    m.add_function(wrap_pyfunction!(read_vcf, m)?)?;
    Ok(())
}
