use std::io::{BufReader, BufRead};
use std::path::{Path, PathBuf};
use std::fs::File;
use anyhow::Result;

use crate::vcf::records::{VcfRecord, VcfRecords};
use crate::vcf::metadata::InfoMeta;
use crate::vcf::elements::*;

use super::metadata::FormatMeta;

pub struct InMemoryVcfRecords {
    records: Vec<VcfRecord>,
}

impl VcfRecords for InMemoryVcfRecords {
    fn from_path<P: AsRef<Path>>(path: &P) -> Result<Self> {
        let mut records: Vec<VcfRecord> = Vec::new();
        let mut info_metas: Vec<InfoMeta> = Vec::new();
        let mut format_metas: Vec<FormatMeta> = Vec::new(); 
        let mut sample_ids: Vec<String> = Vec::new();
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
            let pos = row_values[1].parse::<i128>()?;
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
            records.push(
                VcfRecord::new(chrom, pos, id, ref_, alt, qual, filter, info, format)
            )

        }

        Ok(InMemoryVcfRecords {
            records,
        })
    }

    fn get_by_id(&self, id: &str) -> Option<VcfRecord> {
        for record in &self.records {
            if record.id() == id {
                return Some(record.clone())
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_path() {
        let records = InMemoryVcfRecords::from_path(&"examples/demo_merge/test.merge.manta.vcf".to_string()).unwrap();
        let second = records.records.get(1).unwrap();
        assert_eq!(second.id(), "MD1");
        assert_eq!(second.chrom().get_chr_pretended_value(), "chr1");
        assert_eq!(second.alt()[0], AltElement::AngleBracket("DEL".to_string()));
        assert_eq!(second.info().get("SVTYPE").unwrap()[0], VcfData::String("DEL".to_string()));
        assert_eq!(second.format().get("mouse01_N").unwrap().get("SR").unwrap()[0], VcfData::Integer(10));
    }
}
