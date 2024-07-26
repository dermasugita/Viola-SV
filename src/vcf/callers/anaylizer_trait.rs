use crate::vcf::{elements::{AltElement, End, VcfData}, records::VcfRecord};
use crate::vcf::callers::tables::*;


pub trait RecordAnalyzer {
    fn positions(&self, record: &VcfRecord) -> PositionsTableRow;
    fn infos(&self, record: &VcfRecord) -> (StrInfoTable, IntInfoTable, FloatInfoTable, FlagInfoTable);
    fn formats(&self, record: &VcfRecord) -> (StrFormatTable, IntFormatTable, FloatFormatTable);
}
