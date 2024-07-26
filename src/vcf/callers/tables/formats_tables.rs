use polars::prelude::*;
use pyo3::prelude::*;
use anyhow::{Result, anyhow};

#[pyclass]
#[derive(Clone)]
pub struct FormatTables {
    #[pyo3(get)]
    pub str_format: StrFormatTable,
    #[pyo3(get)]
    pub int_format: IntFormatTable,
    #[pyo3(get)]
    pub float_format: FloatFormatTable,
}
impl FormatTables {
    pub fn new() -> Self {
        Self {
            str_format: StrFormatTable::new(),
            int_format: IntFormatTable::new(),
            float_format: FloatFormatTable::new(),
        }
    }
    pub fn from_tables(str_format: StrFormatTable, int_format: IntFormatTable, float_format: FloatFormatTable) -> Self {
        Self {
            str_format,
            int_format,
            float_format,
        }
    }
    pub fn push(&mut self, formats: (StrFormatTable, IntFormatTable, FloatFormatTable)) {
        self.str_format = self.str_format.cat(&formats.0).expect("Error concatenating StrFormatTables");
        self.int_format = self.int_format.cat(&formats.1).expect("Error concatenating IntFormatTables");
        self.float_format = self.float_format.cat(&formats.2).expect("Error concatenating FloatFormatTables");
    }
    pub fn cat(&self, other: &Self) -> Result<Self> {
        let str_table = self.str_format.cat(&other.str_format)?;
        let int_table = self.int_format.cat(&other.int_format)?;
        let float_table = self.float_format.cat(&other.float_format)?;
        Ok(Self::from_tables(str_table, int_table, float_table))
    }
}
pub struct StrFormatTableRow {
    pub id: String,
    pub sample: String,
    pub format: String,
    pub value_idx: u32,
    pub value: String,
}

#[pyclass]
#[derive(Clone)]
pub struct StrFormatTable {
    _df: DataFrame,
}
impl StrFormatTable {
    pub fn new() -> Self {
        let _df = df!(
            "id" => Vec::<String>::new(),
            "sample" => Vec::<String>::new(),
            "format" => Vec::<String>::new(),
            "value_idx" => Vec::<u32>::new(),
            "value" => Vec::<String>::new()
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn from_vec(vec: &[StrFormatTableRow]) -> Self {
        let mut id = Vec::new();
        let mut sample = Vec::new();
        let mut format = Vec::new();
        let mut value_idx = Vec::new();
        let mut value = Vec::new();
        for row in vec {
            id.push(row.id.clone());
            sample.push(row.sample.clone());
            format.push(row.format.clone());
            value_idx.push(row.value_idx);
            value.push(row.value.clone());
        }
        let _df = df!(
            "id" => id,
            "sample" => sample,
            "format" => format,
            "value_idx" => value_idx,
            "value" => value
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn push(&mut self, row: StrFormatTableRow) {
        let push_df = Self::from_vec(&[row]);
        self._df.vstack_mut(&push_df._df).expect("Error pushing row");
    }
    pub fn get_by_ids(&self, ids: &[&str]) -> Result<Self> {
        let mask = self._df.column("id")?.equal(&Series::new("id", ids))?;
        Ok(
            Self {
                _df: self._df.filter(&mask)?
            }
        )
    }
    pub fn cat(&self, other: &Self) -> Result<Self> {
        match self._df.vstack(other.get_df()) {
            Ok(df) => Ok(Self { _df: df }),
            Err(e) => Err(anyhow!("Error concatenating DataFrames: {}", e)),

        }
    }
    pub fn get_df(&self) -> &DataFrame {
        &self._df
    }
    pub fn get_df_mut(&mut self) -> &mut DataFrame {
        &mut self._df
    }
}
#[pymethods]
impl StrFormatTable {
    pub fn get_columns(&self) -> Vec<String> {
        self._df.get_column_names().into_iter().map(|s| s.to_string()).collect()
    }
    #[getter]
    pub fn id(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("id").expect("id not found")
            .str().expect("id value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("ID not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn sample(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("sample").expect("sample not found")
            .str().expect("sample value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Sample not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn format(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("format").expect("format not found")
            .str().expect("format value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Format not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn value_idx(&self) -> PyResult<Vec<u32>> {
        let ret = Vec::from(self._df.column("value_idx").expect("value_idx not found")
            .u32().expect("value_idx value cannot be converted to u32"))
            .into_iter()
            .map(|s| s.expect("Value index not found")).collect();
        Ok(ret)
    }
    #[getter]
    pub fn value(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("value").expect("value not found")
            .str().expect("value value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Value not found").to_string()).collect();
        Ok(ret)
    }
}

pub struct IntFormatTableRow {
    pub id: String,
    pub sample: String,
    pub format: String,
    pub value_idx: u32,
    pub value: i64,
}
#[pyclass]
#[derive(Clone)]
pub struct IntFormatTable {
    _df: DataFrame,
}

impl IntFormatTable {
    pub fn new() -> Self {
        let _df = df!(
            "id" => Vec::<String>::new(),
            "sample" => Vec::<String>::new(),
            "format" => Vec::<String>::new(),
            "value_idx" => Vec::<u32>::new(),
            "value" => Vec::<i64>::new()
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn from_vec(vec: &[IntFormatTableRow]) -> Self {
        let mut id = Vec::new();
        let mut sample = Vec::new();
        let mut format = Vec::new();
        let mut value_idx = Vec::new();
        let mut value = Vec::new();
        for row in vec {
            id.push(row.id.clone());
            sample.push(row.sample.clone());
            format.push(row.format.clone());
            value_idx.push(row.value_idx);
            value.push(row.value);
        }
        let _df = df!(
            "id" => id,
            "sample" => sample,
            "format" => format,
            "value_idx" => value_idx,
            "value" => value
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn push(&mut self, row: IntFormatTableRow) {
        let push_df = Self::from_vec(&[row]);
        self._df.vstack_mut(&push_df._df).expect("Error pushing row");
    }
    pub fn get_by_ids(&self, ids: &[&str]) -> Result<Self> {
        let mask = self._df.column("id")?.equal(&Series::new("id", ids))?;
        Ok(
            Self {
                _df: self._df.filter(&mask)?
            }
        )
    }
    pub fn cat(&self, other: &Self) -> Result<Self> {
        match self._df.vstack(other.get_df()) {
            Ok(df) => Ok(Self { _df: df }),
            Err(e) => Err(anyhow!("Error concatenating DataFrames: {}", e)),

        }
    }
    pub fn get_df(&self) -> &DataFrame {
        &self._df
    }
    pub fn get_df_mut(&mut self) -> &mut DataFrame {
        &mut self._df
    }
}
#[pymethods]
impl IntFormatTable {
    pub fn get_columns(&self) -> Vec<String> {
        self._df.get_column_names().into_iter().map(|s| s.to_string()).collect()
    }
    #[getter]
    pub fn id(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("id").expect("id not found")
            .str().expect("id value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("ID not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn sample(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("sample").expect("sample not found")
            .str().expect("sample value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Sample not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn format(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("format").expect("format not found")
            .str().expect("format value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Format not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn value_idx(&self) -> PyResult<Vec<u32>> {
        let ret = Vec::from(self._df.column("value_idx").expect("value_idx not found")
            .u32().expect("value_idx value cannot be converted to u32"))
            .into_iter()
            .map(|s| s.expect("Value index not found")).collect();
        Ok(ret)
    }
    #[getter]
    pub fn value(&self) -> PyResult<Vec<i64>> {
        let ret = Vec::from(self._df.column("value").expect("value not found")
            .i64().expect("value value cannot be converted to i64"))
            .into_iter()
            .map(|s| s.expect("Value not found")).collect();
        Ok(ret)
    }
}

pub struct FloatFormatTableRow {
    pub id: String,
    pub sample: String,
    pub format: String,
    pub value_idx: u32,
    pub value: Option<f32>,
}

#[pyclass]
#[derive(Clone)]
pub struct FloatFormatTable {
    _df: DataFrame,
}
impl FloatFormatTable {
    pub fn new() -> Self {
        let _df = df!(
            "id" => Vec::<String>::new(),
            "sample" => Vec::<String>::new(),
            "format" => Vec::<String>::new(),
            "value_idx" => Vec::<u32>::new(),
            "value" => Vec::<Option<f32>>::new()
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn from_vec(vec: &[FloatFormatTableRow]) -> Self {
        let mut id = Vec::new();
        let mut sample = Vec::new();
        let mut format = Vec::new();
        let mut value_idx = Vec::new();
        let mut value = Vec::new();
        for row in vec {
            id.push(row.id.clone());
            sample.push(row.sample.clone());
            format.push(row.format.clone());
            value_idx.push(row.value_idx);
            value.push(row.value);
        }
        let _df = df!(
            "id" => id,
            "sample" => sample,
            "format" => format,
            "value_idx" => value_idx,
            "value" => value
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn push(&mut self, row: FloatFormatTableRow) {
        let push_df = Self::from_vec(&[row]);
        self._df.vstack_mut(&push_df._df).expect("Error pushing row");
    }
    pub fn get_by_ids(&self, ids: &[&str]) -> Result<Self> {
        let mask = self._df.column("id")?.equal(&Series::new("id", ids))?;
        Ok(
            Self {
                _df: self._df.filter(&mask)?
            }
        )
    }
    pub fn cat(&self, other: &Self) -> Result<Self> {
        match self._df.vstack(other.get_df()) {
            Ok(df) => Ok(Self { _df: df }),
            Err(e) => Err(anyhow!("Error concatenating DataFrames: {}", e)),

        }
    }
    pub fn get_df(&self) -> &DataFrame {
        &self._df
    }
    pub fn get_df_mut(&mut self) -> &mut DataFrame {
        &mut self._df
    }
}

#[pymethods]
impl FloatFormatTable {
    pub fn get_columns(&self) -> Vec<String> {
        self._df.get_column_names().into_iter().map(|s| s.to_string()).collect()
    }
    #[getter]
    pub fn id(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("id").expect("id not found")
            .str().expect("id value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("ID not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn sample(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("sample").expect("sample not found")
            .str().expect("sample value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Sample not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn format(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("format").expect("format not found")
            .str().expect("format value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Format not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn value_idx(&self) -> PyResult<Vec<u32>> {
        let ret = Vec::from(self._df.column("value_idx").expect("value_idx not found")
            .u32().expect("value_idx value cannot be converted to u32"))
            .into_iter()
            .map(|s| s.expect("Value index not found")).collect();
        Ok(ret)
    }
    #[getter]
    pub fn value(&self) -> PyResult<Vec<Option<f64>>> {
        let res = Vec::from(self._df.column("value").expect("value not found").f64().expect("value value cannot be converted to f64"));
        Ok(res)
    }
}