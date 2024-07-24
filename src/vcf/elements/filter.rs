#[derive(Debug, Clone)]
pub struct Filter {
    value: Option<String>
}

impl Filter {
    pub fn from_str(value: &str) -> Self {
        if value == "." {
            Self {
                value: None
            }
        } else {
            Self {
                value: Some(value.to_string())
            }
        }
    }
}