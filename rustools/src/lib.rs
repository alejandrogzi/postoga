use pyo3::prelude::*;

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::string::String;

use flate2::write::GzEncoder;
use flate2::Compression;
use natord::compare;
use rayon::prelude::*;

pub mod bed;
pub mod codon;
pub mod lines;
pub mod utils;

#[pyfunction]
#[pyo3(signature = (bed, output, isoforms=".", gz = true, no_gene = false))]
fn convert(
    py: Python,
    bed: PyObject,
    output: PyObject,
    isoforms: Option<&str>,
    gz: bool,
    no_gene: bool,
) -> PyResult<PathBuf> {
    let bed = bed.extract::<PathBuf>(py)?;
    let isoforms = isoforms.map(PathBuf::from);
    let output = output.extract::<PathBuf>(py)?;

    convert_to_gxf(&bed, isoforms, &output, gz, no_gene);
    Ok(output)
}

#[pymodule]
#[pyo3(name = "rustools")]
fn rustools(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(convert, m)?)?;
    m.add_function(wrap_pyfunction!(extract_seqs, m)?)?;
    Ok(())
}

fn convert_to_gxf(
    bed: &PathBuf,
    isoforms: Option<PathBuf>,
    output: &PathBuf,
    gz: bool,
    no_gene: bool,
) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_cpus::get() as usize)
        .build()
        .unwrap();

    let filetype = output.extension().unwrap().to_str().unwrap();
    let imap = if !no_gene {
        utils::get_isoforms(&isoforms.unwrap())
    } else {
        HashMap::new()
    };

    let bed = utils::bed_reader(&bed);
    let gene_track = utils::custom_par_parse(&bed).unwrap_or_else(|_| {
        let message = format!("Error parsing bed file");
        panic!("{}", message);
    });

    let results = bed
        .par_iter()
        .filter_map(|record| match filetype {
            "gff" | "gff3" => to_gxf::<b'='>(record, &imap).ok(),
            "gtf" => to_gxf::<b' '>(record, &imap).ok(),
            _ => None,
        })
        .flatten()
        .collect::<Vec<_>>();

    let mut blocks = match filetype {
        "gff" | "gff3" => utils::combine_maps_par::<b'='>(&imap, &gene_track),
        "gtf" => utils::combine_maps_par::<b' '>(&imap, &gene_track),
        "gz" => {
            match output
                .file_stem()
                .and_then(|s| s.to_str())
                .and_then(|s| s.rsplit('.').next())
            {
                Some("gff") | Some("gff3") => utils::combine_maps_par::<b'='>(&imap, &gene_track),
                Some("gtf") => utils::combine_maps_par::<b' '>(&imap, &gene_track),
                _ => panic!("ERROR: Invalid output file type -> {:?}. Tried to match its extension and failed", output),
            }
        }
        _ => panic!("Invalid output file type"),
    };
    blocks.extend(results);

    blocks.par_sort_unstable_by(|a, b| {
        let chr_cmp = compare(&a.0, &b.0);
        if chr_cmp == std::cmp::Ordering::Equal {
            a.2.cmp(&b.2)
        } else {
            chr_cmp
        }
    });

    let writer_boxed: Box<dyn Write> = if gz {
        let file = File::create(&output).unwrap();
        let encoder = GzEncoder::new(file, Compression::default());
        Box::new(BufWriter::new(encoder))
    } else {
        let file = File::create(&output).unwrap();
        Box::new(BufWriter::new(file))
    };

    let mut writer = writer_boxed;
    for entry in &blocks {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t{}",
            entry.0, "bed2gxf", entry.1, entry.2, entry.3, entry.4, entry.5, entry.6
        )
        .unwrap();
    }
}

#[inline(always)]
fn to_gxf<const SEP: u8>(
    bedline: &bed::BedRecord,
    isoforms: &HashMap<String, String>,
) -> Result<Vec<(String, String, u32, u32, String, String, String)>, Box<dyn Error>> {
    let mut result: Vec<(String, String, u32, u32, String, String, String)> = Vec::new();

    let gene = if !isoforms.is_empty() {
        match isoforms.get(&bedline.name) {
            Some(g) => g,
            None => {
                panic!("Gene {} not found in isoforms file.", bedline.name);
            }
        }
    } else {
        &bedline.name
    };

    let fcodon = codon::first_codon(bedline)
        .unwrap_or_else(|| panic!("No start codon found for {}.", bedline.name));
    let lcodon = codon::last_codon(bedline).unwrap_or_else(|| {
        panic!("No stop codon found for {}.", bedline.name);
    });

    let frames = bedline.get_frames();
    let cds_end: u32 = if bedline.strand == "+" && codon::codon_complete(&lcodon) {
        move_pos(bedline, lcodon.end, -3)
    } else {
        bedline.cds_end
    };

    let cds_start = if bedline.strand == "-" && codon::codon_complete(&fcodon) {
        move_pos(bedline, fcodon.start, 3)
    } else {
        bedline.cds_start
    };

    match SEP {
        b'=' => lines::build_gff_line(
            bedline,
            gene,
            "transcript",
            bedline.tx_start,
            bedline.tx_end,
            3,
            -1,
            &mut result,
        ),
        b' ' => lines::build_gtf_line(
            bedline,
            gene,
            "transcript",
            bedline.tx_start,
            bedline.tx_end,
            3,
            -1,
            &mut result,
        ),
        _ => panic!("ERROR: Invalid format. Could not convert BED record to GTF line"),
    }

    for i in 0..bedline.exon_count as usize {
        match SEP {
            b'=' => lines::build_gff_line(
                bedline,
                gene,
                "exon",
                bedline.exon_start[i],
                bedline.exon_end[i],
                3,
                i as i16,
                &mut result,
            ),
            b' ' => lines::build_gtf_line(
                bedline,
                gene,
                "exon",
                bedline.exon_start[i],
                bedline.exon_end[i],
                3,
                i as i16,
                &mut result,
            ),
            _ => panic!("ERROR: Invalid format. Could not convert BED record to GTF line"),
        }

        if cds_start < cds_end {
            lines::write_features::<SEP>(
                i,
                bedline,
                gene,
                cds_start,
                cds_end,
                frames[i] as u32,
                &mut result,
            );
        }
    }

    if bedline.strand != "-" {
        if codon::codon_complete(&fcodon) {
            lines::write_codon::<SEP>(bedline, gene, "start_codon", fcodon, &mut result);
        }
        if codon::codon_complete(&lcodon) {
            lines::write_codon::<SEP>(bedline, gene, "stop_codon", lcodon, &mut result);
        }
    } else {
        if codon::codon_complete(&lcodon) {
            lines::write_codon::<SEP>(bedline, gene, "start_codon", lcodon, &mut result);
        }
        if codon::codon_complete(&fcodon) {
            lines::write_codon::<SEP>(bedline, gene, "stop_codon", fcodon, &mut result);
        }
    }

    Ok(result)
}

fn move_pos(record: &bed::BedRecord, pos: u32, dist: i32) -> u32 {
    let mut pos = pos;
    assert!(record.tx_start <= pos && pos <= record.tx_end);

    let mut exon_index = record
        .exon_start
        .iter()
        .zip(record.exon_end.iter())
        .position(|(start, end)| pos >= *start && pos <= *end)
        .unwrap_or_else(|| {
            let message = format!("Position {} not in exons.", pos);
            panic!("{}", message);
        }) as i16;

    let mut steps = dist.abs();
    let direction = if dist >= 0 { 1 } else { -1 };

    while steps > 0 {
        let (exon_start, exon_end) = (
            record.exon_start[exon_index as usize],
            record.exon_end[exon_index as usize],
        );

        if pos >= exon_start && pos <= exon_end {
            pos += direction as u32;
            steps -= 1;
        } else if direction >= 0 {
            exon_index += 1;
            if (exon_index as usize) < record.exon_count as usize {
                pos = record.exon_start[exon_index as usize];
            }
        } else {
            exon_index -= 1;
            if exon_index >= 0 {
                pos = record.exon_end[exon_index as usize] - 1;
                steps -= 1;
            }
        }
    }
    if steps > 0 {
        panic!("can't move {} by {}", pos, dist);
    }
    pos
}

#[pyfunction]
#[pyo3(signature = (bed, fasta, hint="query", output=None))]
fn extract_seqs(
    py: Python,
    bed: PyObject,
    fasta: PyObject,
    hint: &str,
    output: Option<PyObject>,
) -> PyResult<PathBuf> {
    let bed = bed.extract::<PathBuf>(py)?;
    let fasta = fasta.extract::<PathBuf>(py)?;
    let output = match output {
        Some(output) => output.extract::<PathBuf>(py)?,
        None => {
            let mut output = fasta.clone();
            output.set_extension("filtered.fa");
            output
        }
    };
    let hint = utils::Hint::from_str(hint);

    let records = utils::reader(&bed)?;
    let txs = utils::extract_tx_from_bed(records.as_str());

    let seqs = utils::reader(&fasta)?;
    let fa = utils::build_fasta_hash(seqs.as_bytes(), hint).unwrap_or_else(|_| {
        panic!(
            "ERROR: Could not build FASTA hash from {}.",
            fasta.display()
        );
    });

    let mut out = BufWriter::new(File::create(&output)?);

    for tx in txs {
        let seq = fa.get(tx).unwrap_or_else(|| {
            panic!("ERROR: {} not found in FASTA file.", tx);
        });
        writeln!(out, ">{}", tx)?;
        writeln!(out, "{}", seq)?;
    }

    Ok(output)
}
