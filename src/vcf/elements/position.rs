use derive_getters::Getters;

use crate::vcf::elements::contig::Contig;

#[derive(Getters, Clone, Debug, PartialEq, Eq)]
pub struct Position {
    chrom: Contig,
    pos: i128,
}

impl Position {
    pub fn new(chrom: Contig, pos: i128) -> Self {
        Self {
            chrom,
            pos,
        }
    }
    pub fn from_colon_separated(s: &str) -> Self {
        let mut iter = s.split(':');
        let chrom = Contig::new(iter.next().unwrap());
        let pos = iter.next().unwrap().parse().unwrap();
        Self {
            chrom,
            pos,
        }
    }
    pub fn get_colon_separated(&self) -> String {
        format!("{}:{}", self.chrom.get_chr_pretended_value(), self.pos)
    }
}