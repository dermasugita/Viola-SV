use std::collections::HashMap;
use anyhow::{Result, anyhow};
use crate::vcf::metadata::{VcfKV, VcfNestedKV};
pub enum ABCVcfHeader {
    KV(VcfKV),
    NestedKV(VcfNestedKV),
}

#[derive(Eq, PartialEq, Hash, Copy, Clone, Debug)]
enum State {
    Start,
    ReadKey1,
    PatternDivision,
    ReadKey2,
    ReadValue,
    InString,
    End
}

struct HeaderParserFSM {
    state: State,
    key1: String,
    key2: String,
    value: String,
    kv: HashMap<String, String>
}

impl HeaderParserFSM {
    fn new() -> Self {
        Self {
            state: State::Start,
            key1: String::new(),
            key2: String::new(),
            value: String::new(),
            kv: HashMap::new()
        }
    }

    fn run(&mut self, s: &str) -> Result<ABCVcfHeader> {
        for c in s.chars() {
            match self.state {
                State::Start => {
                    match c {
                        '#' => {},
                        _ => {
                            self.key1.push(c);
                            self.state = State::ReadKey1;
                        }
                    }
                },
                State::ReadKey1 => {
                    match c {
                        '=' => {
                            self.state = State::PatternDivision;
                        },
                        _ => {
                            self.key1.push(c);
                        }
                    }
                },
                State::PatternDivision => {
                    match c {
                        '<' => {
                            self.state = State::ReadKey2;
                        },
                        _ => {
                            self.value.push(c);
                            self.state = State::ReadValue;
                        }
                    }
                }
                State::ReadKey2 => {
                    match c {
                        '=' => {
                            self.state = State::ReadValue;
                        },
                        _ => {
                            self.key2.push(c);
                        }
                    }
                },
                State::ReadValue => {
                    match c {
                        ',' => {
                            self.kv.insert(self.key2.clone(), self.value.clone());
                            self.key2.clear();
                            self.value.clear();
                            self.state = State::ReadKey2;
                        },
                        '>' => {
                            self.kv.insert(self.key2.clone(), self.value.clone());
                            self.key2.clear();
                            self.value.clear();
                            self.state = State::End;
                        },
                        '"' => {
                            self.state = State::InString;
                        },
                        _ => {
                            self.value.push(c);
                        }
                    }
                },
                State::InString => {
                    match c {
                        '"' => {
                            self.state = State::ReadValue;
                        },
                        _ => {
                            self.value.push(c);
                        }
                    }
                },
                State::End => {
                    match c {
                        ' ' => { },
                        _ => { return Err(anyhow!("Invalid character found after the end of the header: {}", c)); }
                    }
                }
            }
        }

        if self.kv.len() > 0 {
            return Ok(ABCVcfHeader::NestedKV(VcfNestedKV::new(self.key1.clone(), self.kv.clone())));
        }
        Ok(ABCVcfHeader::KV(VcfKV::new(self.key1.clone(), self.value.clone())))
        
    }
}

pub fn parse_nested_vcf_header(s: &str) -> Result<HashMap<String, String>> {
    let mut fsm = HeaderParserFSM::new();
    match fsm.run(s) {
        Ok(ABCVcfHeader::NestedKV(kv)) => {
            Ok(kv.get_kv())
        },
        Ok(ABCVcfHeader::KV(_)) => {
            Err(anyhow!("The header is not nested."))
        },
        Err(e) => {
            Err(e)
        }
    }
}
pub fn parse_vcf_header(s: &str) -> Result<ABCVcfHeader> {
    let mut fsm = HeaderParserFSM::new();
    fsm.run(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_header_parser_fsm() {
        let mut fsm = HeaderParserFSM::new();
        let s = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
        let result = fsm.run(s).unwrap();
        let result = match result {
            ABCVcfHeader::NestedKV(nkv) => nkv.get_kv(),
            _ => panic!("Invalid result.")
        };
        assert_eq!(result.get("ID").unwrap(), "GT");
        assert_eq!(result.get("Number").unwrap(), "1");
        assert_eq!(result.get("Type").unwrap(), "String");
        assert_eq!(result.get("Description").unwrap(), "Genotype");
    }
    
    #[test]
    fn test_header_parser_with_comma() {
        let mut fsm = HeaderParserFSM::new();
        let s = "##INFO=<ID=JUNCTION_SOMATICSCORE,Number=1,Type=Integer,Description=\"If the SV junctino is part of an EVENT (ie. a multi-adjacency variant), this field provides the SOMATICSCORE value for the adjacency in question only\">";
        let result = fsm.run(s).unwrap();
        let result = match result {
            ABCVcfHeader::NestedKV(nkv) => nkv.get_kv(),
            _ => panic!("Invalid result.")
        };
        assert_eq!(result.get("ID").unwrap(), "JUNCTION_SOMATICSCORE");
        assert_eq!(result.get("Number").unwrap(), "1");
        assert_eq!(result.get("Type").unwrap(), "Integer");
        assert_eq!(result.get("Description").unwrap(), "If the SV junctino is part of an EVENT (ie. a multi-adjacency variant), this field provides the SOMATICSCORE value for the adjacency in question only");
    }

    #[test]
    fn test_header_parser_kv() {
        let mut fsm = HeaderParserFSM::new();
        let s = "##fileformat=VCFv4.1";
        let result = fsm.run(s).unwrap();
        let result = match result {
            ABCVcfHeader::KV(kv) => kv.get_value(),
            _ => panic!("Invalid result.")
        };
        assert_eq!(result, "VCFv4.1");
    }
}