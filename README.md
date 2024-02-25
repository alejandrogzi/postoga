# postoga

The post-TOGA processing pipeline.

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![version](https://img.shields.io/badge/version-0.8.0--devel-orange)

<p align="center">
    <img width=700 align="center" src="./supply/postoga_logo_git.png" >
</p>

<!-- <img src="./supply/postoga_report.png" align="center"/> -->


## Usage

To use postoga, just do:

```bash
# clone the repository
git clone --recursive https://github.com/alejandrogzi/postoga.git
cd postoga

# create an environment (recommended)
conda env create -f env.yml
conda activate postoga

# OR call configure.sh to install dependencies
./configure.sh

# run a test to confirm functionality
./test.sh
```

If you see something like this at then end, postoga is ready!:

```text
####################################
postoga: the post-TOGA processing pipeline
version: 0.5.0-devel
commit: 67ecb4e
branch: master

[XXXXXX] - INFO: postoga started!
[XXXXXX] - INFO: running in mode base with arguments: {'mode': 'base', 'path': './supply/test/', 'by_class': None, 'by_rel': None, 'threshold': '0.5', 'to': 'gff', 'assembly_qual': './supply/Ancestral_placental.txt'}
[XXXXXX] - INFO: found 10 projections, 10 unique transcripts, 10 unique genes
```

> [!NOTE]
>
> postoga is dependent from TOGA. Any changes in TOGA will have a repercusion here. If you found any bug/errors, please report them here.
> This project is in constant development, any desired features are welcome!

Here is a descrption of postoga features:

```text
usage: postoga.py [-h] {base,haplotype}

positional arguments:
  {base,haplotype}  Select mode
    base            Base mode
    haplotype       Haplotype mode

postoga.py base [-h] -p PATH [-bc BY_CLASS] [-br BY_REL] [-th THRESHOLD] [-par PARALOG] [-to {gtf,gff,bed}] [-aq ASSEMBLY_QUAL] [-sp {human,mouse,chicken}]

optional arguments:
  -h, --help            Show this help message and exit
  -p PATH, --path PATH  Path to TOGA results directory (all postoga results will be send to that location)
  -bc BY_CLASS, --by-class BY_CLASS
                        Filter parameter to only include certain orthology classes (I, PI, UL, M, PM, L, UL), specified with commas (e.g. I,PI,UL).
  -br BY_REL, --by-rel BY_REL
                        Filter parameter to only include certain orthology relationships (o2o, o2m, m2m, m2m, o2z), specified with commas (e.g. o2o,o2m).
  -th THRESHOLD, --threshold THRESHOLD
                        Filter parameter to preserve orthology scores greater or equal to a given threshold (0.0 - 1.0)
  -to {gtf,gff,bed}, --to {gtf,gff,bed}
                        Specify the conversion format for .bed (query_annotation/filtered) file (gtf, gff3) or just keep it as .bed (bed)
  -aq ASSEMBLY_QUAL, --assembly_qual ASSEMBLY_QUAL
                        Calculate assembly quality based on a list of genes provided by the user (default: Ancestral_placental.txt)
  -sp {human,mouse,chicken}, --species {human,mouse,chicken}
                        Species name to be used as a reference for the assembly quality calculation (default: human)
  -src {ensembl,gene_name,entrez}, --source {ensembl,gene_name,entrez}
                        Source of the ancestral gene names (default: ENSG)
  -phy {mammals,birds}, --phylo {mammals,birds}
                        Phylogenetic group of your species (default: mammals)
  -s, --skip            Skip steps 2, 3, and 4 and only filter the .bed file
  -par PARALOG, --paralog PARALOG
                        Filter parameter to preserve transcripts with paralog projection probabilities less or equal to a given threshold (0.0 - 1.0)

postoga.py haplotype [-h] -hp HAPLOTYPE_PATH [-r RULE] [-s {query,loss}]

optional arguments:
  -h, --help            Display help message
  -hp HAPLOTYPE_PATH, --haplotype_path HAPLOTYPE_PATH
                        Path to TOGA results directories separated by commas (path1,path2,path3)
  -r RULE, --rule RULE  Rule to merge haplotype assemblies (default: I>PI>UL>L>M>PM>PG>NF)
  -s {query,loss}, --source {query,loss}
                        Source of the haplotype classes (default: loss)
```


## What's new on version 0.8.0-devel

- Adds `--paralog`, a new argument to preserve transcripts with paralog projections probabilities less or equal to a given threshold.
- Implements `--skip` to only use postoga as a filtering tool, skipping additional post-processing steps.
- Adds `bed` as an option to `--to` to avoid the conversion to gtf/gff.
- Modifies `logger.py` to send logging information to stdout (this is also written to `postoga.log`)

