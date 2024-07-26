#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SupportedCallers {
    Manta,
    Delly,
    Lumpy,
    Gridss,
}

impl SupportedCallers {
    pub fn from_str(s: &str) -> Self {
        match s {
            "manta" => SupportedCallers::Manta,
            "delly" => SupportedCallers::Delly,
            "lumpy" => SupportedCallers::Lumpy,
            "gridss" => SupportedCallers::Gridss,
            _ => panic!("Invalid caller name: {}", s),
        }
    }
}