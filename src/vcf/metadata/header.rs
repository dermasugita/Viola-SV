
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