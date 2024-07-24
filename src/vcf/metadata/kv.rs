
use std::collections::HashMap;
use derive_getters::Getters;

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