use std::char;

use crate::vcf::records::*;
use crate::vcf::elements::*;
use crate::vcf::callers::*;
pub struct MantaAnalyzer;
impl RecordAnalyzer for MantaAnalyzer {
    fn positions(&self, record: &VcfRecord) -> PositionsTableRow {
        let id = record.id().to_string();
        let chrom1 = record.chrom().value().clone();
        let mut pos1 = *record.pos();
        let mut chrom2 = None;
        let mut pos2 = None;
        let mut strand1 = '+';
        let mut strand2 = None;
        let qual = *record.qual();
        let ref_ = record.ref_().to_string();
        let mut svtype = "".to_string();

        let alt = record.alt();


        match &alt[0] {
            AltElement::Breakend(alt) => {
                svtype = "BND".to_string();
                if let End::ThreePrime = alt.this_end() { strand1 = '+'; } else { strand1 = '-'; }
                if let Some(position) = alt.other_end_pos() {
                    chrom2 = Some(position.chrom().value().clone());
                    pos2 = Some(*position.pos());
                    if let Some(End::ThreePrime) = alt.other_end_extend_to() { strand2 = Some('-'); } else { strand2 = Some('+'); }
                }
            },
            AltElement::AngleBracket(alt) => {
                let end = record.info().get("END").expect(&format!("END not found on record: id {}", &record.id()));
                let end = if let VcfData::Integer(end) = end[0] {
                    end as i64
                } else {
                    panic!("END is not an integer on record id: {}", &record.id());
                };
                match alt.as_str() {
                    "DEL" => {
                        chrom2 = Some(chrom1.clone());
                        pos2 = Some(end + 1);
                        strand1 = '+';
                        strand2 = Some('-');
                        svtype = "DEL".to_string();
                    },
                    "DUP" | "DUP:TANDEM" => {
                        chrom2 = Some(chrom1.clone());
                        pos2 = Some(end + 1);
                        strand1 = '-';
                        strand2 = Some('+');
                        svtype = "DUP".to_string();
                    },
                    "INV" => {
                        let inv3 = record.info().get("INV3");
                        chrom2 = Some(chrom1.clone());
                        svtype = "INV".to_string();
                        if let None = inv3 {
                            pos2 = Some(end);
                            strand1 = '+';
                            strand2 = Some('+');
                        } else {
                            pos1 = pos1 + 1;
                            pos2 = Some(end + 1);
                            strand1 = '-';
                            strand2 = Some('-');
                        }
                    },
                    _ => {
                        panic!("Unsupported ALT: {}", alt);
                    }
                }
            },
            _ => {
                panic!("Unsupported ALT: {:?}", alt);
            }
        };

        PositionsTableRow {
            id,
            chrom1,
            pos1,
            chrom2,
            pos2,
            strand1,
            strand2,
            qual,
            ref_,
            alt: alt[0].to_string(),
            svtype,
        }
    }
    fn infos(&self, record: &VcfRecord) -> (StrInfoTable, IntInfoTable, FloatInfoTable, FlagInfoTable) {
        let mut char_info = Vec::new();
        let mut string_info = Vec::new();
        let mut int_info = Vec::new();
        let mut float_info = Vec::new();
        let mut flag_info = Vec::new();
        for (info_id, value_vec) in record.info().get_hashmap() {
            for (idx, value) in value_vec.iter().enumerate() {
                match value {
                    VcfData::Character(value) => {
                        char_info.push(StrInfoTableRow {
                            id: record.id().to_string(),
                            key: info_id.to_string(),
                            value_idx: idx as u32,
                            value: value.to_string(),
                        });
                    },
                    VcfData::String(value) => {
                        string_info.push(StrInfoTableRow {
                            id: record.id().to_string(),
                            key: info_id.to_string(),
                            value_idx: idx as u32,
                            value: value.to_string(),
                        });
                    },
                    VcfData::Integer(value) => {
                        int_info.push(IntInfoTableRow {
                            id: record.id().to_string(),
                            key: info_id.to_string(),
                            value_idx: idx as u32,
                            value: *value as i64,
                        });
                    },
                    VcfData::Float(value) => {
                        float_info.push(FloatInfoTableRow {
                            id: record.id().to_string(),
                            key: info_id.to_string(),
                            value_idx: idx as u32,
                            value: *value,
                        });
                    },
                    VcfData::Flag(_) => {
                        flag_info.push(FlagInfoTableRow {
                            id: record.id().to_string(),
                            key: info_id.to_string(),
                            value_idx: idx as u32,
                        });
                    },
                }
            }
        }
        (
            StrInfoTable::from_vec(&string_info),
            IntInfoTable::from_vec(&int_info),
            FloatInfoTable::from_vec(&float_info),
            FlagInfoTable::from_vec(&flag_info)
        )
    }
    fn formats(&self, record: &VcfRecord) -> (StrFormatTable, IntFormatTable, FloatFormatTable) {
        let mut str_format = Vec::new();
        let mut int_format = Vec::new();
        let mut float_format = Vec::new();
        for (sample_id, format_unit) in record.format().get_hashmap() {
            for (format_id, format_values) in format_unit.get_hashmap() {
                for (value_idx, format_values) in format_values.into_iter().enumerate() {
                    match format_values {
                        VcfData::Character(value) => {
                            str_format.push(StrFormatTableRow {
                                id: record.id().clone(),
                                sample: sample_id.clone(),
                                format: format_id.to_string(),
                                value_idx: value_idx as u32,
                                value: value.to_string(),
                            })
                        },
                        VcfData::String(value) => {
                            str_format.push(StrFormatTableRow {
                                id: record.id().clone(),
                                sample: sample_id.clone(),
                                format: format_id.to_string(),
                                value_idx: value_idx as u32,
                                value: value.to_string(),
                            })
                        },
                        VcfData::Integer(value) => {
                            int_format.push(IntFormatTableRow {
                                id: record.id().clone(),
                                sample: sample_id.clone(),
                                format: format_id.to_string(),
                                value_idx: value_idx as u32,
                                value: *value as i64,
                            })
                        },
                        VcfData::Float(value) => {
                            float_format.push(FloatFormatTableRow {
                                id: record.id().clone(),
                                sample: sample_id.clone(),
                                format: format_id.to_string(),
                                value_idx: value_idx as u32,
                                value: *value,
                            })
                        },
                        VcfData::Flag(_) => {
                            panic!("Flag format is not specified in VCF specification. record_id: {}, format: {}", record.id(), format_id.to_string());
                        }
                    }
                } 
            }
        }
        (
            StrFormatTable::from_vec(&str_format),
            IntFormatTable::from_vec(&int_format),
            FloatFormatTable::from_vec(&float_format),
        )
    }
}