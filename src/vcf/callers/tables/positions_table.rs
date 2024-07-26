use pyo3::prelude::*;
use anyhow::{Result, anyhow};
use indexmap::IndexSet;
use polars::prelude::*;

#[derive(Clone)]
pub struct PositionsTableRow {
    pub id: String,
    pub chrom1: String,
    pub pos1: i64,
    pub chrom2: Option<String>,
    pub pos2: Option<i64>,
    pub strand1: char,
    pub strand2: Option<char>,
    pub qual: Option<f32>,
    pub ref_: String,
    pub alt: String,
    pub svtype: String,
}
#[pyclass]
#[derive(Clone)]
pub struct PositionsTable {
    _df: DataFrame,
}

impl PositionsTable {
    pub fn new() -> Self {
        let _df = df!(
            "id" => Vec::<String>::new(),
            "chrom1" => Vec::<String>::new(),
            "pos1" => Vec::<i64>::new(),
            "chrom2" => Vec::<Option<String>>::new(),
            "pos2" => Vec::<Option<i64>>::new(),
            "strand1" => Vec::<String>::new(),
            "strand2" => Vec::<Option<String>>::new(),
            "qual" => Vec::<f32>::new(),
            "ref_" => Vec::<String>::new(),
            "alt" => Vec::<String>::new(),
            "svtype" => Vec::<String>::new(),
        ).expect("Error creating DataFrame");
        Self {
            _df
        }
    }
    pub fn from_vec(vec: &[PositionsTableRow]) -> Self {
        let mut _id = IndexSet::new();
        let mut chrom1 = Vec::new();
        let mut pos1 = Vec::new();
        let mut chrom2 = Vec::new();
        let mut pos2 = Vec::new();
        let mut strand1 = Vec::new();
        let mut strand2 = Vec::new();
        let mut qual = Vec::new();
        let mut ref_ = Vec::new();
        let mut alt = Vec::new();
        let mut svtype = Vec::new();
        for row in vec {
            _id.insert(row.id.clone());
            chrom1.push(row.chrom1.clone());
            pos1.push(row.pos1);
            chrom2.push(row.chrom2.clone());
            pos2.push(row.pos2);
            strand1.push(row.strand1.to_string());
            strand2.push(match row.strand2{
                Some(s) => Some(s.to_string()),
                None => None,
            });
            qual.push(row.qual);
            ref_.push(row.ref_.clone());
            alt.push(row.alt.clone());
            svtype.push(row.svtype.clone());
        }
        let _df = df!(
            "id" => _id.into_iter().collect::<Vec<String>>(),
            "chrom1" => chrom1,
            "pos1" => pos1,
            "chrom2" => chrom2,
            "pos2" => pos2,
            "strand1" => strand1,
            "strand2" => strand2,
            "qual" => qual,
            "ref_" => ref_,
            "alt" => alt,
            "svtype" => svtype,
        ).expect("Error creating DataFrame");
        Self {
            _df
        }
    }
    pub fn push(&mut self, row: PositionsTableRow) {
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
    pub fn get_df(&self) -> &DataFrame {
        &self._df
    }
    pub fn get_df_mut(&mut self) -> &mut DataFrame {
        &mut self._df
    }
}
#[pymethods]
impl PositionsTable {
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
    pub fn chrom1(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("chrom1").expect("chrom1 not found")
            .str().expect("chrom1 value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Chrom1 not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn pos1(&self) -> PyResult<Vec<i64>> {
        let ret = Vec::from(self._df.column("pos1").expect("pos1 not found")
            .i64().expect("pos1 value cannot be converted to i64"))
            .into_iter()
            .map(|s| s.expect("Pos1 not found")).collect();
        Ok(ret)
    }
    #[getter]
    pub fn chrom2(&self) -> PyResult<Vec<Option<String>>> {
        let ret = Vec::from(self._df.column("chrom2").expect("chrom2 not found")
            .str().expect("chrom2 value cannot be converted to string"))
            .into_iter()
            .map(|s| match s {
                Some(s) => Some(s.to_string()),
                None => None,
            }).collect();
        Ok(ret)
    }
    #[getter]
    pub fn pos2(&self) -> PyResult<Vec<Option<i64>>> {
        let ret = Vec::from(self._df.column("pos2").expect("pos2 not found")
            .i64().expect("pos2 value cannot be converted to i64"))
            .into_iter()
            .map(|s| match s {
                Some(s) => Some(s),
                None => None,
            }).collect();
        Ok(ret)
    }
    #[getter]
    pub fn strand1(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("strand1").expect("strand1 not found")
            .str().expect("strand1 value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Strand1 not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn strand2(&self) -> PyResult<Vec<Option<String>>> {
        let ret = Vec::from(self._df.column("strand2").expect("strand2 not found")
            .str().expect("strand2 value cannot be converted to string"))
            .into_iter()
            .map(|s| match s {
                Some(s) => Some(s.to_string()),
                None => None,
            }).collect();
        Ok(ret)
    }
    #[getter]
    pub fn qual(&self) -> PyResult<Vec<Option<f32>>> {
        let ret = Vec::from(self._df.column("qual").expect("qual not found")
            .f32().expect("qual value cannot be converted to f32"))
            .into_iter()
            .map(|s| match s {
                Some(s) => Some(s),
                None => None,
            }).collect();
        Ok(ret)
    }
    #[getter]
    pub fn ref_(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("ref_").expect("ref_ not found")
            .str().expect("ref_ value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Ref_ not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn alt(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("alt").expect("alt not found")
            .str().expect("alt value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Alt not found").to_string()).collect();
        Ok(ret)
    }
    #[getter]
    pub fn svtype(&self) -> PyResult<Vec<String>> {
        let ret = Vec::from(self._df.column("svtype").expect("svtype not found")
            .str().expect("svtype value cannot be converted to string"))
            .into_iter()
            .map(|s| s.expect("Svtype not found").to_string()).collect();
        Ok(ret)
    }
}