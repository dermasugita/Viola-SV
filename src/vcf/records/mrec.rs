use std::io::{BufReader, BufRead};
use std::path::{Path, PathBuf};
use std::fs::File;
use indexmap::IndexMap;
use anyhow::{Result, anyhow};

use crate::vcf::callers::*;
use crate::vcf::records::{VcfRecord, VcfRecords};
use crate::vcf::metadata::*;
use crate::vcf::elements::*;


pub struct InMemoryVcfRecords {
    records: IndexMap<String, VcfRecord>,
    current_index: usize,
    positions: PositionsTable,
    infos: InfoTables,
    formats: FormatTables
}
impl InMemoryVcfRecords {
    fn _compute_positions_table<Analyzer: RecordAnalyzer>(&self, analyzer: &Analyzer) -> PositionsTable {
        let mut positions_table = PositionsTable::new();
        for record in self.records.values() {
            let positions = analyzer.positions(record);
            positions_table.push(positions);
        }
        positions_table
    }
}

impl VcfRecords for InMemoryVcfRecords {
    fn from_path<P: AsRef<Path>>(path: &P, caller: &SupportedCallers) -> Result<Self> {
        let mut records: IndexMap<String, VcfRecord> = IndexMap::new();
        let mut info_metas: Vec<InfoMeta> = Vec::new();
        let mut format_metas: Vec<FormatMeta> = Vec::new(); 
        let mut sample_ids: Vec<String> = Vec::new();
        let mut positions: PositionsTable = PositionsTable::new();
        let mut infos: InfoTables = InfoTables::new();
        let mut formats: FormatTables = FormatTables::new();
        for line in BufReader::new(File::open(path)?).lines() {
            let line = line?;
            if line.starts_with("##INFO") {
                let info_meta = InfoMeta::from_str(&line)?;
                info_metas.push(info_meta);
                continue;
            }
            if line.starts_with("##FORMAT") {
                let format_meta = FormatMeta::from_str(&line)?;
                format_metas.push(format_meta);
                continue;
            }
            if line.starts_with("#CHROM") {
                let headers = line.split("\t").collect::<Vec<&str>>();
                if headers.len() >= 9 {
                    sample_ids = headers[9..].iter().map(|s| s.to_string()).collect::<Vec<String>>();
                }
                continue;
            }
            if line.starts_with("#") {
                continue;
            }
            let row_values = line.split("\t").collect::<Vec<&str>>();
            let chrom = Contig::new(row_values[0]);
            let pos = row_values[1].parse::<i64>()?;
            let id = row_values[2].to_string();
            let ref_ = DNABaseSeq::from_str(row_values[3])?;
            let alt = row_values[4].split(",").map(|s| AltElement::from_str(s).unwrap()).collect::<Vec<AltElement>>();
            let qual = if row_values[5] == "." { None } else { Some(row_values[5].parse::<f32>()?) };
            let filter = row_values[6].split(',').map(|s| Filter::from_str(s)).collect::<Vec<Filter>>();
            let info = Info::from_str(row_values[7], &info_metas)?;
            let mut format = Format::new();
            let format_guide = row_values[8];
            for (sample_id, format_values) in row_values[9..].iter().zip(sample_ids.iter()) {
                let format_unit = FormatUnit::from_str(format_guide, sample_id, &format_metas)?;
                format.insert(format_values.to_string(), format_unit);
            }
            let record = VcfRecord::new(chrom, pos, id.clone(), ref_, alt, qual, filter, info, format);
            records.insert(
                id,
                record.clone());
            match caller {
                SupportedCallers::Manta => {
                    let analyzer = MantaAnalyzer;
                    let positions_row = analyzer.positions(&record);
                    let infos_row = analyzer.infos(&record);
                    let formats_row = analyzer.formats(&record);
                    positions.push(positions_row);
                    infos.push(infos_row);
                    formats.push(formats_row);
                },
                _ => {return Err(anyhow!("Unsupported caller"));}
            }
        }

        Ok(InMemoryVcfRecords {
            records,
            current_index: 0,
            positions,
            infos,
            formats
        })
    }

    fn get_by_ids(&self, ids: &[&str]) -> Result<Vec<VcfRecord>> {
        let mut records = Vec::new();
        let mut not_found = Vec::new();
        let mut has_not_found = false;
        for id in ids {
            if has_not_found {
                if let None = self.records.get(*id) {not_found.push(id.to_string()); continue;}
            }
            match self.records.get(*id) {
                Some(record) => records.push(record.clone()),
                None => {not_found.push(id.to_string()); has_not_found = true;}
            }
        }
        if not_found.len() > 0 {
            return Err(anyhow!("The following ids are not found: {:?}", not_found));
        }
        Ok(records)
    }
    fn remove_by_id(&mut self, id: &str) -> Option<VcfRecord> {
        self.records.swap_remove(id)
    }
    fn count(&self) -> usize {
        self.records.len()
    }
    fn get_positions_table(&self) -> &PositionsTable {
        &self.positions
    }
    fn get_info_tables(&self) -> &InfoTables {
        &self.infos
    }
    fn get_format_tables(&self) -> &FormatTables {
        &self.formats
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_path() {
        let records = InMemoryVcfRecords::from_path(&"examples/demo_merge/test.merge.manta.vcf".to_string(), &SupportedCallers::Manta).unwrap();
        let second = records.records.get_index(1).unwrap().1;
        assert_eq!(second.id(), "MD1");
        assert_eq!(second.chrom().get_chr_pretended_value(), "chr1");
        assert_eq!(second.alt()[0], AltElement::AngleBracket("DEL".to_string()));
        assert_eq!(second.info().get("SVTYPE").unwrap()[0], VcfData::String("DEL".to_string()));
        assert_eq!(second.format().get("mouse01_N").unwrap().get("SR").unwrap()[0], VcfData::Integer(10));
        assert_eq!(second.ref_().to_string(), "N");
    }
}
impl Iterator for InMemoryVcfRecords {
    type Item = VcfRecord;
    fn next(&mut self) -> Option<Self::Item> {
        match self.records.get_index(self.current_index) {
            Some(record) => {
                self.current_index += 1;
                Some(record.1.clone())
            },
            None => None
        }
    }

}

