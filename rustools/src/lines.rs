use crate::bed::BedRecord;
use crate::codon::*;
use std::cmp::{max, min};
use std::fmt::Write;

pub fn build_gtf_line(
    record: &BedRecord,
    gene: &String,
    gene_type: &str,
    exon_start: u32,
    exon_end: u32,
    frame: u32,
    exon: i16,
    result: &mut Vec<(String, String, u32, u32, String, String, String)>,
) {
    assert!(record.tx_start < record.tx_end);

    let phase = match frame {
        0 => "0",
        1 => "2",
        2 => "1",
        _ => ".",
    };

    let mut attr = format!("gene_id \"{}\"; transcript_id \"{}\";", gene, record.name);

    if exon >= 0 {
        let (exon_id, nexon) = if record.strand == "+" {
            let exon_id = exon + 1;
            (exon_id as u16, exon + 1)
        } else {
            let exon_id = record.exon_count - exon as u16;
            (exon_id, exon_id as i16)
        };

        write!(
            attr,
            " exon_number \"{}\"; exon_id \"{}.{}\";",
            nexon, record.name, exon_id
        )
        .expect("Failed to write exon information");
    }

    result.push((
        record.chrom.clone(),
        gene_type.to_string(),
        exon_start + 1,
        exon_end,
        record.strand.clone(),
        phase.to_string(),
        attr,
    ));
}

pub fn build_gff_line(
    record: &BedRecord,
    gene: &String,
    gene_type: &str,
    exon_start: u32,
    exon_end: u32,
    frame: u32,
    exon: i16,
    result: &mut Vec<(String, String, u32, u32, String, String, String)>,
) {
    assert!(record.tx_start < record.tx_end);

    let phase = match frame {
        0 => "0",
        1 => "2",
        2 => "1",
        _ => ".",
    };

    let mut attr = String::new();

    if gene_type == "transcript" {
        attr.push_str(&format!(
            "ID={};Parent={};gene_id={};transcript_id={}",
            record.name, gene, gene, record.name
        ));
    } else {
        if exon >= 0 {
            let (exon_id, nexon) = if record.strand == "+" {
                let exon_id = exon + 1;
                (exon_id as u16, exon + 1)
            } else {
                let exon_id = record.exon_count - exon as u16;
                (exon_id, exon_id as i16)
            };

            attr.push_str(&format!(
                "ID={}:{}.{};Parent={};gene_id={};transcript_id={};exon_number={}",
                gene_type, record.name, exon_id, record.name, gene, record.name, nexon
            ));
        } else {
            let prefix = match gene_type {
                "five_prime_utr" => "5UTR",
                "three_prime_utr" => "3UTR",
                _ => panic!("Invalid gene type"),
            };

            attr.push_str(&format!(
                "ID={}:{};Parent={};gene_id={};transcript_id={}",
                prefix, record.name, record.name, gene, record.name
            ));
        }
    }

    result.push((
        record.chrom.clone(),
        gene_type.to_string(),
        exon_start + 1,
        exon_end,
        record.strand.clone(),
        phase.to_string(),
        attr,
    ));
}

pub fn write_features<const SEP: u8>(
    i: usize,
    record: &BedRecord,
    gene: &String,
    cds_start: u32,
    cds_end: u32,
    frame: u32,
    result: &mut Vec<(String, String, u32, u32, String, String, String)>,
) {
    let exon_start = record.exon_start[i];
    let exon_end = record.exon_end[i];

    if record.cds_start < exon_end && exon_start < record.cds_end {
        let start = max(exon_start, cds_start);
        let end = min(exon_end, cds_end);
        if start < end {
            match SEP {
                b'=' => build_gff_line(record, gene, "CDS", start, end, frame, i as i16, result),
                b' ' => build_gtf_line(record, gene, "CDS", start, end, frame, i as i16, result),
                _ => panic!("Invalid separator"),
            }
        }
    }
}

pub fn write_codon<const SEP: u8>(
    record: &BedRecord,
    gene: &String,
    gene_type: &str,
    codon: Codon,
    result: &mut Vec<(String, String, u32, u32, String, String, String)>,
) {
    match SEP {
        b'=' => build_gff_line(
            record,
            gene,
            gene_type,
            codon.start,
            codon.end,
            0,
            codon.index as i16,
            result,
        ),
        b' ' => build_gtf_line(
            record,
            gene,
            gene_type,
            codon.start,
            codon.end,
            0,
            codon.index as i16,
            result,
        ),
        _ => panic!("Invalid separator"),
    }

    if codon.start2 < codon.end2 {
        match SEP {
            b'=' => build_gff_line(
                record,
                gene,
                gene_type,
                codon.start,
                codon.end,
                codon.start2,
                (codon.end - codon.start) as i16,
                result,
            ),
            b' ' => build_gtf_line(
                record,
                gene,
                gene_type,
                codon.start,
                codon.end,
                codon.start2,
                (codon.end - codon.start) as i16,
                result,
            ),
            _ => panic!("Invalid separator"),
        }
    }
}
