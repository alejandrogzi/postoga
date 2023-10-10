# postoga

The post-TOGA processing pipeline.

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![version](https://img.shields.io/badge/version-0.6.0--devel-orange)

<img src="./supply/postoga_logo_git.png" width=500/>
![Report](https://github.com/alejandrogzi/postoga/blob/master/supply/postoga_report.png)

## Usage

To use postoga, just do:

```bash
# clone the repository
git clone https://github.com/alejandrogzi/postoga.git
cd postoga

# call configure.sh to install dependencies
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

Here is a descrption of postoga features:

```text
usage: postoga.py [-h] {base,haplotype}

positional arguments:
  {base,haplotype}  Select mode
    base            Base mode
    haplotype       Haplotype mode

postoga.py base [-h] -p PATH [-bc BY_CLASS] [-br BY_REL] [-th THRESHOLD] -to {gtf,gff} [-aq ASSEMBLY_QUAL] [-sp {human,mouse,chicken}]

optional arguments:
  -h, --help            Display help message
  -p PATH, --path PATH  Path to TOGA results directory
  -bc BY_CLASS, --by-class BY_CLASS
                        Filter parameter to only include certain orthology classes (I, PI, UL, M, PM, L, UL)
  -br BY_REL, --by-rel BY_REL
                        Filter parameter to only include certain orthology relationships (o2o, o2m, m2m, m2m, o2z)
  -th THRESHOLD, --threshold THRESHOLD
                        Filter parameter to preserve orthology scores greater or equal to a given threshold (0.0 - 1.0)
  -to {gtf,gff}, --to {gtf,gff}
                        Specify the conversion format for .bed (query_annotation/filtered) file (gtf, gff3)
  -aq ASSEMBLY_QUAL, --assembly_qual ASSEMBLY_QUAL
                        Calculate assembly quality based on a list of genes provided by the user (default: Ancestral_placental.txt)
  -sp {human,mouse,chicken}, --species {human,mouse,chicken}
                        Species name to be used as a reference for the assembly quality calculation (default: human)

postoga.py haplotype [-h] -hp HAPLOTYPE_PATH [-r RULE] [-s {query,loss}]

optional arguments:
  -h, --help            Display help message
  -hp HAPLOTYPE_PATH, --haplotype_path HAPLOTYPE_PATH
                        Path to TOGA results directories separated by commas (path1,path2,path3)
  -r RULE, --rule RULE  Rule to merge haplotype assemblies (default: I>PI>UL>L>M>PM>PG>NF)
  -s {query,loss}, --source {query,loss}
                        Source of the haplotype classes (default: loss)
```


## What's new on version 0.6.0-devel



- Fixed unnecessary imports and comments in all modules.
- Added `plotter.py`, the plotter module of postoga.
- postoga now reports findings automatically and save them under `POSTOGA_REPORT.pdf`
- Modules have been updated to synchronize with plotter module.
- Added plotter-dependent project-wide constants
- Implemented `get_stats_from_bed` under `filter_query_annotation`, to quickly extract query stats 



## TO DO's

- [x] Build logger/version/constants scripts

- [x] postoga STEP 1: automate conversion from bed to gtf/gff -> bed2gtf will implement a sorting leaving dependecy on gtfsort (something already implemented in bed2gff

- [x] test model with different naming nomenclatures (most recent TOGA versions)

- [ ] Handle possible warnings of data (low number of intact genes, etc) within log file (?)

- [x] Check/install/test requeriments (install_dependencies.py) -> rust, cargo, binaries (bed2gtf, bed2gff), compleasm (?)

- [x] postoga STEP 2: assembly stats re-implementation -> https://github.com/hillerlab/TOGA/blob/master/supply/TOGA_assemblyStats.py

- [x] Build Graph class automating plotting reports -> need to define canvas structure and plot types to a default report (relevant variables)

- [ ] postoga STEP 3: alignment stats -> rust binding + ortholog lengths (?)

- [ ] postoga STEP 4: implement [compleasm](https://github.com/huangnengCSU/compleasm) an efficient alternative to BUSCO -> could be compared with assembly_stats.py -> then compared with an existente database of mammalian/placental genomes (?)

- [ ] ask Bogdan/Michael/Scott feedback and ideas