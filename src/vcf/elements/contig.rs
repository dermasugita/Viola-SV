use derive_getters::Getters;
#[derive(Getters, Clone, Debug, PartialEq, Eq)]
pub struct Contig {
    value: String,
    chr_pretended: bool
}
impl Contig {
    pub fn new(value: &str) -> Self {
        let mut chr_pretended = false;
        if value.starts_with("chr") {
            chr_pretended = true;
        }
        Self {
            value: value.to_string(),
            chr_pretended
        }
    }
    pub fn get_chr_pretended_value(&self) -> String {
        if self.chr_pretended {
            self.value.clone()
        } else {
            format!("chr{}", self.value)
        }
    }
    pub fn get_chr_removed_value(&self) -> String {
        if self.chr_pretended {
            self.value.trim_start_matches("chr").to_string()
        } else {
            self.value.clone()
        }
    }
}