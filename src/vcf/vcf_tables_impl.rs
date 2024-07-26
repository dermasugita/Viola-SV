use anyhow::Result;
use crate::vcf::vcf::Vcf;
use crate::vcf::records::VcfRecords;
use crate::vcf::callers::*;

impl<T: VcfRecords> Vcf<T> {
    pub fn get_positions_table(&mut self) -> &PositionsTable {
        self.get_records().get_positions_table()
    }
    pub fn get_info_tables(&mut self)
        -> &InfoTables {
        self.get_records().get_info_tables()
    }
    pub fn get_format_tables(&mut self)
        -> &FormatTables {
        self.get_records().get_format_tables()
    }
    
}