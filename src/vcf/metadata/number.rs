use anyhow::{Result, anyhow};

#[derive(Copy, Clone)]
pub enum VcfNumber {
    Number(i32),
    A,
    G,
    R,
    Dot,
}
impl VcfNumber {
    pub fn from_str(s: &str) -> Result<Self> {
        match s {
            "A" => Ok(VcfNumber::A),
            "G" => Ok(VcfNumber::G),
            "R" => Ok(VcfNumber::R),
            "." => Ok(VcfNumber::Dot),
            _ => {
                match s.parse::<i32>() {
                    Ok(n) => Ok(VcfNumber::Number(n)),
                    Err(_) => Err(anyhow!("Unknown VcfNumber has been found. {}", s)),
                }
            }
        }
    }
}