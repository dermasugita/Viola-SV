use pyo3::prelude::*;
use polars::prelude::*;
use anyhow::{Result, anyhow};

#[pyclass]
#[derive(Clone)]
pub struct InfoTables {
    #[pyo3(get)]
    pub str_info: StrInfoTable,
    #[pyo3(get)]
    pub int_info: IntInfoTable,
    #[pyo3(get)]
    pub float_info: FloatInfoTable,
    #[pyo3(get)]
    pub flag_info: FlagInfoTable,
}
impl InfoTables {
    pub fn new() -> Self {
        Self {
            str_info: StrInfoTable::new(),
            int_info: IntInfoTable::new(),
            float_info: FloatInfoTable::new(),
            flag_info: FlagInfoTable::new(),
        }
    }
    pub fn from_tables(str_info: StrInfoTable, int_info: IntInfoTable, float_info: FloatInfoTable, flag_info: FlagInfoTable) -> Self {
        Self {
            str_info,
            int_info,
            float_info,
            flag_info,
        }
    }
    pub fn push(&mut self, infos: (StrInfoTable, IntInfoTable, FloatInfoTable, FlagInfoTable)) {
        self.str_info = self.str_info.cat(&infos.0).expect("Error concatenating StrInfoTables");
        self.int_info = self.int_info.cat(&infos.1).expect("Error concatenating IntInfoTables");
        self.float_info = self.float_info.cat(&infos.2).expect("Error concatenating FloatInfoTables");
        self.flag_info = self.flag_info.cat(&infos.3).expect("Error concatenating FlagInfoTables");
    }
    pub fn cat(&self, other: &Self) -> Result<Self> {
        Ok(
            Self {
                str_info: self.str_info.cat(&other.str_info)?,
                int_info: self.int_info.cat(&other.int_info)?,
                float_info: self.float_info.cat(&other.float_info)?,
                flag_info: self.flag_info.cat(&other.flag_info)?,
            }
        )
    }
}

pub struct StrInfoTableRow {
    pub id: String,
    pub key: String,
    pub value_idx: u32,
    pub value: String,
}

#[pyclass]
#[derive(Clone)]
pub struct StrInfoTable {
    _df: DataFrame,
}
impl StrInfoTable {
    pub fn new() -> Self {
        let _df = df!(
            "id" => Vec::<String>::new(),
            "key" => Vec::<String>::new(),
            "value_idx" => Vec::<u32>::new(),
            "value" => Vec::<String>::new(),
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn from_vec(vec: &[StrInfoTableRow]) -> Self {
        let mut id = Vec::new();
        let mut key = Vec::new();
        let mut value_idx = Vec::new();
        let mut value = Vec::new();
        for row in vec {
            id.push(row.id.clone());
            key.push(row.key.clone());
            value_idx.push(row.value_idx as u32);
            value.push(row.value.to_string());
        }
        let _df = df!(
            "id" => id,
            "key" => key,
            "value_idx" => value_idx,
            "value" => value,
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn push(&mut self, row: StrInfoTableRow) {
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
impl StrInfoTable {
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
    pub fn key(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("key").expect("key not found")
            .str().expect("key value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Key not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn value_idx(&self) -> PyResult<Vec<u32>> {
        let ret = Vec::from(self._df.column("value_idx").expect("value_idx not found")
            .u32().expect("value_idx value cannot be converted to u32"))
            .into_iter()
            .map(|s| s.expect("Value_idx not found")).collect();
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

pub struct IntInfoTableRow {
    pub id: String,
    pub key: String,
    pub value_idx: u32,
    pub value: i64,
}
#[pyclass]
#[derive(Clone)]
pub struct IntInfoTable {
    _df: DataFrame
}
impl IntInfoTable {
    pub fn new() -> Self {
        let _df = df!(
            "id" => Vec::<String>::new(),
            "key" => Vec::<String>::new(),
            "value_idx" => Vec::<u32>::new(),
            "value" => Vec::<i64>::new(),
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn from_vec(vec: &[IntInfoTableRow]) -> Self {
        let mut id = Vec::new();
        let mut key = Vec::new();
        let mut value_idx = Vec::new();
        let mut value = Vec::new();
        for row in vec {
            id.push(row.id.clone());
            key.push(row.key.clone());
            value_idx.push(row.value_idx);
            value.push(row.value);
        }
        let _df = df!(
            "id" => id,
            "key" => key,
            "value_idx" => value_idx,
            "value" => value,
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn push(&mut self, row: IntInfoTableRow) {
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
impl IntInfoTable {
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
    pub fn key(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("key").expect("key not found")
            .str().expect("key value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Key not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn value_idx(&self) -> PyResult<Vec<u32>> {
        let ret = Vec::from(self._df.column("value_idx").expect("value_idx not found")
            .u32().expect("value_idx value cannot be converted to u32"))
            .into_iter()
            .map(|s| s.expect("Value_idx not found")).collect();
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

pub struct FloatInfoTableRow {
    pub id: String,
    pub key: String,
    pub value_idx: u32,
    pub value: Option<f32>,
}
#[pyclass]
#[derive(Clone)]
pub struct FloatInfoTable {
    _df: DataFrame
}
impl FloatInfoTable {
    pub fn new() -> Self {
        let _df = df!(
            "id" => Vec::<String>::new(),
            "key" => Vec::<String>::new(),
            "value_idx" => Vec::<u32>::new(),
            "value" => Vec::<Option<f32>>::new(),
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn from_vec(vec: &[FloatInfoTableRow]) -> Self {
        let mut id = Vec::new();
        let mut key = Vec::new();
        let mut value_idx = Vec::new();
        let mut value = Vec::new();
        for row in vec {
            id.push(row.id.clone());
            key.push(row.key.clone());
            value_idx.push(row.value_idx);
            value.push(row.value);
        }
        let _df = df!(
            "id" => id,
            "key" => key,
            "value_idx" => value_idx,
            "value" => value,
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn push(&mut self, row: FloatInfoTableRow) {
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
impl FloatInfoTable {
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
    pub fn key(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("key").expect("key not found")
            .str().expect("key value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Key not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn value_idx(&self) -> PyResult<Vec<u32>> {
        let ret = Vec::from(self._df.column("value_idx").expect("value_idx not found")
            .u32().expect("value_idx value cannot be converted to u32"))
            .into_iter()
            .map(|s| s.expect("Value_idx not found")).collect();
        Ok(ret)
    }
    #[getter]
    pub fn value(&self) -> PyResult<Vec<Option<f32>>> {
        let ret = Vec::from(self._df.column("value").expect("value not found")
            .f32().expect("value value cannot be converted to f32"));
        Ok(ret)
    }
}


pub struct FlagInfoTableRow {
    pub id: String,
    pub key: String,
    pub value_idx: u32,
}
#[pyclass]
#[derive(Clone)]
pub struct FlagInfoTable {
    _df: DataFrame
}
impl FlagInfoTable {
    pub fn new() -> Self {
        let _df = df!(
            "id" => Vec::<String>::new(),
            "key" => Vec::<String>::new(),
            "value_idx" => Vec::<u32>::new(),
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn from_vec(vec: &[FlagInfoTableRow]) -> Self {
        let mut id = Vec::new();
        let mut key = Vec::new();
        let mut value_idx = Vec::new();
        for row in vec {
            id.push(row.id.clone());
            key.push(row.key.clone());
            value_idx.push(row.value_idx);
        }
        let _df = df!(
            "id" => id,
            "key" => key,
            "value_idx" => value_idx,
        ).expect("Error creating DataFrame");
        Self { _df }
    }
    pub fn push(&mut self, row: FlagInfoTableRow) {
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
impl FlagInfoTable {
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
    pub fn key(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("key").expect("key not found")
            .str().expect("key value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Key not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn value_idx(&self) -> PyResult<Vec<usize>> {
        let ret = Vec::from(self._df.column("value_idx").expect("value_idx not found")
            .u32().expect("value_idx value cannot be converted to u32"))
            .into_iter()
            .map(|s| s.expect("Value_idx not found") as usize).collect();
        Ok(ret)
    }
}