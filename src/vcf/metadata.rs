use std::collections::HashMap;
use derive_getters::Getters;
use url::Url;

#[derive(Copy, Clone)]
pub enum VcfDataType {
    Integer,
    Float,
    Flag,
    Character,
    String,
}

#[derive(Eq, PartialEq, Hash, Clone, Debug)]
pub enum ViolaVcfInfoId {
    CIPOS,
    CIEND,
    IMPRECISE,
    PRECISE,
    SVTYPE,
    Others(String),
}

#[derive(Eq, PartialEq, Hash, Clone, Debug)]
pub enum ViolaVcfFormatId {
    PR,
    SR,
    Others(String),
}

#[derive(Copy, Clone)]
pub enum VcfNumber {
    Number(i32),
    A,
    G,
    R,
    Dot,
}

#[derive(Getters, Clone)]
pub struct VcfKV<K = String, V = String>
where K: Clone, V: Clone {
    key: K,
    value: V,
}
impl<K, V> VcfKV<K, V>
where K: Clone, V: Clone {
    pub fn new(key: K, value: V) -> Self {
        VcfKV {
            key,
            value,
        }
    }
    pub fn get_value(&self) -> V {
        self.value.clone()
    }
}

#[derive(Getters, Clone)]
pub struct VcfNestedKV<K1 = String, K2 = String, V = String> 
where K2: Clone, V: Clone {
    key: K1,
    value: HashMap<K2, V>,
}
impl<K1, K2, V> VcfNestedKV<K1, K2, V>
where K2: Clone, V: Clone {
    pub fn new(key: K1, value: HashMap<K2, V>) -> Self {
        VcfNestedKV {
            key,
            value,
        }
    }
    pub fn get_kv(&self) -> HashMap<K2, V> {
        self.value.clone()
    }
}

#[derive(Getters, Clone)]
pub struct ContigMeta {
    id: String,
    length: Option<i128>,
    assembly: Option<String>,
    url: Option<Url>,
    md5: Option<String>,
    species: Option<String>,
    taxonomy: Option<String>,
}
impl ContigMeta {
    pub fn new(id: String, length: Option<i128>, assembly: Option<String>, url: Option<Url>, md5: Option<String>, species: Option<String>, taxonomy: Option<String>) -> Self {
        ContigMeta {
            id,
            length,
            assembly,
            url,
            md5,
            species,
            taxonomy,
        }
    }

}


#[derive(Getters, Clone)]
pub struct InfoMeta {
    id: ViolaVcfInfoId,
    number: VcfNumber,
    dtype: VcfDataType,
    description: String,
    source: Option<String>,
    version: Option<String>,
}
impl InfoMeta {
    pub fn new(id: ViolaVcfInfoId, number: VcfNumber, dtype: VcfDataType, description: String, source: Option<String>, version: Option<String>) -> Self {
        InfoMeta {
            id,
            number,
            dtype,
            description,
            source,
            version,
        }
    }
}

#[derive(Getters, Clone)]
pub struct FilterMeta {
    id: String,
    description: String,
}
impl FilterMeta {
    pub fn new(id: String, description: String) -> Self {
        FilterMeta {
            id,
            description,
        }
    }
}

#[derive(Getters, Clone)]
pub struct FormatMeta {
    id: ViolaVcfFormatId,
    number: VcfNumber,
    dtype: VcfDataType,
    description: String,
}

impl FormatMeta {
    pub fn new(id: ViolaVcfFormatId, number: VcfNumber, dtype: VcfDataType, description: String) -> Self {
        FormatMeta {
            id,
            number,
            dtype,
            description,
        }
    }
}

#[derive(Getters, Clone)]
pub struct AltMeta {
    id: String,
    description: String,
}
impl AltMeta {
    pub fn new(id: String, description: String) -> Self {
        AltMeta {
            id,
            description,
        }
    }
}


// Below 3 structs are used for sample meta data
#[derive(Clone, Getters)]
pub struct MetaId {
    name: String,
}
impl MetaId {
    pub fn new(name: String) -> Self {
        MetaId {
            name,
        }
    }
}

#[derive(Clone, Getters)]
pub struct Meta {
    id: MetaId,
    dtype: VcfDataType,
    number: VcfNumber,
    values: Vec<String>,
    description: String,
}
impl Meta {
    pub fn new(id: MetaId, dtype: VcfDataType, number: VcfNumber, values: Vec<String>, description: String) -> Self {
        Meta {
            id,
            dtype,
            number,
            values,
            description,
        }
    }
}

pub struct SampleMeta {
    id: String,
    description: String,
    kv: Vec<VcfKV<MetaId, String>>,
}
impl SampleMeta {
    pub fn new(id: String, description: String, kv: Vec<VcfKV<MetaId, String>>) -> Self {
        SampleMeta {
            id,
            description,
            kv,
        }
    }
}


pub enum VcfColumnRequiredHeader {
    CHROM,
    POS,
    ID,
    REF,
    ALT,
    QUAL,
    FILTER,
    INFO,
}

pub enum VcfColumnOptionalHeader {
    FORMAT,
    SAMPLE(String),
}