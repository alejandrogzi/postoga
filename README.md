> [!WARNING]
>
> postoga is dependent from [TOGA](https://github.com/hillerlab/TOGA). Any changes in TOGA will have a repercusion here. If you found any bug/errors, please report them here.
> This project is in constant development, any desired features are welcome!


# postoga

The post-[TOGA](https://github.com/hillerlab/TOGA) processing pipeline.

![version](https://img.shields.io/badge/version-0.10-orange)

<p align="center">
    <img width=700 align="center" src="./supply/postoga_logo_git.png" >
</p>

<!-- <img src="./supply/postoga_report.png" align="center"/> -->

## What's new on version 0.10

> - Re-implementation of postoga to match TOGA2.0 output
> - Includes self-owned rust bed to gtf coverters through rustools
> - Adds a reusable PostogaLogger plus `-L/--level` to switch between debug, info, warn, or fully silent runs
> - `--depure` now cleans previous POSTOGA_* runs from the parent output directory before starting
> - Introduces `setup_postoga.sh` to set up the `uv` virtualenv, install postoga, and run `maturin develop` automatically
> - Filtering emits discard counts and summary stats in both stdout and `postoga.log`

## Usage

To use postoga, just:

Clone the repository
```bash
# clone the repository
git clone --recursive https://github.com/alejandrogzi/postoga.git
cd postoga
```

Activate the environment and configure binaries
```bash
# Recommended: bootstrap everything (uv venv + pip install + maturin develop)
./setup_postoga.sh

# Manual alternative using Hatch targets
pip install hatch
make configure
```

Here is a descrption of postoga features:

> [!TIP]
>
> If the only thing you want to do is apply some filters to a TOGA result or convert results to GTF/GFF files, I recommend the following command:
>
> ```bash
> ./postoga.py \
> --togadir /your/TOGA/dir \
> --outdir /your/out/dir \
> -bc [YOUR CLASSES] \
> -br [YOUR STATUS] \
> -bs [YOUR THRESHOLD] \
> -bp [YOUR PARALOG SCORE] \
> -to [YOUR FORMAT GTF/GFF] \
> --depure \
> -L info
> ```

```text
usage: postoga.py [-h] --togadir TOGADIR [-bc ORTHOLOGY_STATUS]
                       [-br ORTHOLOGY_CLASS] [-bs ORTHOLOGY_SCORE]
                       [-to {gtf,gff,bed}] [-tg {bed,utr}]
                       [-bp MIN_PARALOG_SCORE] [-w WITH_ISOFORMS]
                       [-ext [{query,reference}]] [--outdir OUTDIR]
                       [--only-table] [--only-convert]
                       [-L {debug,info,warn,off}] [--depure] [--version]

optional arguments:
  -h, --help            show this help message and exit
  --togadir TOGADIR, -td TOGADIR
                        Path to TOGA results directory
  -bc ORTHOLOGY_STATUS, --by-orthology-status ORTHOLOGY_STATUS
                        Filter parameter to only include certain orthology
                        classes (FI, I, PI, UL, M, PM, L, UL)
  -br ORTHOLOGY_CLASS, --by-orthology-class ORTHOLOGY_CLASS
                        Filter parameter to only include certain orthology
                        relationships (o2o, o2m, m2m, m2m, o2z)
  -bs ORTHOLOGY_SCORE, --by-orthology-score ORTHOLOGY_SCORE
                        Filter parameter to preserve orthology scores greater
                        or equal to a given threshold (0.0 - 1.0)
  -to {gtf,gff,bed}, --to {gtf,gff,bed}
                        Specify the conversion format for .bed
                        (query_annotation/filtered) file (gtf, gff3) or just
                        keep it as .bed (bed)
  -tg {bed,utr}, --target {bed,utr}
                        Specify the .bed input file to used by the program
  -bp MIN_PARALOG_SCORE, --by-paralog-score MIN_PARALOG_SCORE
                        Filter parameter to preserve transcripts with paralog
                        projection probabilities less or equal to a given
                        threshold (0.0 - 1.0)
  -w WITH_ISOFORMS, --with-isoforms WITH_ISOFORMS
                        Path to a custom isoform table (default: None)
  -ext [{query,reference}], --extract [{query,reference}]
                        Flag or option to extract sequences (only codon and
                        protein alignments) from the filtered genes. Can be
                        'query', 'reference', or just set as a flag (default:
                        False). When used as a flag extracting 'query'
                        sequences is assumed.
  --outdir OUTDIR, -o OUTDIR
                        Path to posTOGA output directory
  --only-table, -ot    Only produce the toga.table file
  --only-convert, -oc  Only convert the toga.table file to gtf/gff
  -L {debug,info,warn,off}, --level {debug,info,warn,off}
                        Logging verbosity (debug, info, warn, off)
  --depure, -d          Remove any trace of other postoga runs/files
  --version, -v         show this help message and exit
```
