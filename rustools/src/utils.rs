use crate::bed::BedRecord;

use rayon::prelude::*;

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read};
use std::path::PathBuf;

pub fn bed_reader(file: &PathBuf) -> Vec<BedRecord> {
    let bed = reader(file).unwrap();
    let records = parallel_parse(&bed).unwrap();
    records
}

pub fn get_isoforms(f: &PathBuf) -> HashMap<String, String> {
    let file = reader(f).unwrap_or_else(|_| panic!("ERROR: Isoforms file not found. Provide a valid isoform file or run with --no-gene flag."));
    let pairs = parallel_hash_rev(file);

    if pairs.len() == 0 {
        println!(
            "{} {}",
            "ERROR:", "BED file could not be converted. Please check your isoforms file."
        );
        std::process::exit(1);
    }

    pairs
}

pub fn reader(file: &PathBuf) -> io::Result<String> {
    let mut file = File::open(file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    Ok(contents)
}

pub fn parallel_hash<'a>(s: &'a str) -> HashMap<String, String> {
    s.par_lines()
        .filter_map(|line| {
            let mut words = line.split_whitespace();
            if let Some(fw) = words.next() {
                if let Some(sw) = words.next() {
                    return Some((fw.to_owned(), sw.to_owned()));
                }
            }
            // if the line doesn’t have two words, ignore it and return None.
            None
        })
        .collect()
}

pub fn parallel_hash_rev(s: String) -> HashMap<String, String> {
    s.par_lines()
        .filter_map(|line| {
            let mut words = line.split_whitespace();
            if let Some(fw) = words.next() {
                if let Some(sw) = words.next() {
                    return Some((sw.to_owned(), fw.to_owned()));
                }
            }
            None
        })
        .collect()
}

pub fn parallel_parse<'a>(s: &'a str) -> Result<Vec<BedRecord>, &'static str> {
    let records: Result<Vec<BedRecord>, &'static str> =
        s.par_lines().map(|line| BedRecord::parse(line)).collect();

    records
}

pub fn custom_par_parse(
    records: &Vec<BedRecord>,
) -> Result<HashMap<String, (String, u32, u32, String)>, &'static str> {
    let gene_coordinates = records
        .into_par_iter()
        .fold(
            || HashMap::new(),
            |mut acc: HashMap<String, (String, u32, u32, String)>, record| {
                acc.entry(record.name.clone()).or_insert((
                    record.chrom.clone(),
                    record.tx_start,
                    record.tx_end,
                    record.strand.clone(),
                ));
                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut a: HashMap<String, (String, u32, u32, String)>, b| {
                for (key, (chrom, start, end, strand)) in b {
                    a.entry(key).or_insert((chrom, start, end, strand));
                }
                a
            },
        );
    Ok(gene_coordinates)
}

pub fn combine_maps_par<const SEP: u8>(
    isoforms: &HashMap<String, String>,
    gene_track: &HashMap<String, (String, u32, u32, String)>,
) -> Vec<(String, String, u32, u32, String, String, String)> {
    let coords = isoforms
        .par_iter()
        .fold(
            || HashMap::new(),
            |mut acc: HashMap<String, (String, u32, u32, String)>, (transcript, gene)| {
                if let Some(&(ref chrom, start, end, ref strand)) = gene_track.get(transcript) {
                    let entry = acc.entry(gene.clone()).or_insert((
                        chrom.to_string(),
                        start,
                        end,
                        strand.to_string(),
                    ));
                    entry.1 = entry.1.min(start); // update min start
                    entry.2 = entry.2.max(end); // update max end
                }
                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut a, b| {
                for (gene, (chrom, start, end, strand)) in b {
                    let entry = a.entry(gene).or_insert((chrom, start, end, strand));
                    entry.1 = entry.1.min(start); // update min start
                    entry.2 = entry.2.max(end); // update max end
                }
                a
            },
        );

    let lines = coords
        .par_iter()
        .map(|(gene, (chrom, start, end, strand))| {
            let fmt = match SEP {
                b'=' => format!("ID={};gene_id={}", gene, gene),
                b' ' => format!("gene_id \"{}\";", gene),
                _ => format!("gene_id={}", gene),
            };
            (
                chrom.to_string(),
                "gene".to_string(),
                start + 1,
                *end,
                strand.to_string(),
                ".".to_string(),
                fmt,
            )
        })
        .collect();
    lines
}
