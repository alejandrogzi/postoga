## postoga v.0.4.0-devel

- splited main functions in ordered modules following TogaDir class layout (`filter_bed` will be improved in v.0.5.0-devel)
- improved logging (specifies some stats about annotation or filtered annotation [if user choose to filter it]) + connect() method to call it within a module (`bed2gtf`/`bed2gff` log info bothers the structure of postoga log -> will be fixed in v.0.5.0-devel)
- included `transcripts.quality.tsv` in query_table() and calculates some stats about it (included in logging)
- re-implemented `TOGA_assemblyStats.py` to estimate annotation completeness in single mode and haplotype-resolved assemblies -> `--assembly_qual` (to specify genes db -> tab-separated file [.txt,.tsv,.csv,...]; default ancestral); `--haplotype_path` (to specify path of both haplotypes) -> the equivalent of `TOGA_assemblyStats.py -m stats` && `TOGA_assemblyStats.py -m stats -ances ./ances.txt` is completed, `TOGA_assemblyStats.py -m merge` will be available in v.0.5.0-devel.
- postoga now automates: 1) dependencies [v.0.5.0-devel], 2) filtering, 3) conversion, 4) assembly quality estimation (just 1/4 quality steps postoga will have)
- shell is now `utils.py` (postoga base utility module)
- added https://github.com/hillerlab/TOGA/blob/master/TOGAInput/human_hg38/Ancestral_placental.txt to `./supply/Ancestral_placental.txt`

## postoga v.0.5.0-devel

- postoga fully works now on two modes `postoga base` & `postoga haplotype`
- base/haplo branches initiate own logs (haplo branch chooses the first path of the list provided)
- fixed filter_bed() module -> handles different filter combinations (`--by_class`, `--by_rel`, `--threshold`)
- removed unnecessary imports in all modules
- `TOGA_assemblyStats.py -m merge` is now fully functional in postoga under `postoga haplotypes -hpath path1,path2,path3 --source [query, loss] --rule I>PI>UL>L>M>PM>PG>abs`.
- logger.py is automatically updated with the current version
- postoga now automates installing requirements (python/rust) through `./configure sh`
- implemented `test.sh` to make an initial test with random data in `./supply/test`

## postoga v.0.6.0-devel

- Fixed unnecessary imports and comments in all modules.
- Added `plotter.py`, the plotter module of postoga.
- postoga now reports findings automatically and save them under `POSTOGA_REPORT.pdf`
- Modules have been updated to synchronize with plotter module.
- Added plotter-dependent project-wide constants
- Implemented `get_stats_from_bed` under `filter_query_annotation`, to quickly extract query stats

## postoga v.0.7.0-devel

- Now postoga evaluates all the non-overlapping exon lenghts of your annotation and outputs additional stats about it
- Instead of implementing a way to perform BUSCO completeness from scratch, postoga uses a set of curated DBs to evaluate your assembly completeness
- Introducing --source. Now you can choose which gene nomeclature background postoga should use: ensmebl IDs, gene names or Entrez IDs.
- Introducing --phylo. Now you can choose the phylo group of your species to be use a set of BUSCO DBs. See arguments.
- Some changes from the last PR where rolled back. This version of postoga is up to date with TOGA's output dir structure.
- The quality of transcript is no longer considered in this release. postoga now uses orthology prediction scores as an implicit source for quality.

# postoga v.0.8.0-devel

- Adds `--paralog`, a new argument to preserve transcripts with paralog projections probabilities less or equal to a given threshold.
- Implements `--skip` to only use postoga as a filtering tool, skipping additional post-processing steps.
- Adds `bed` as an option to `--to` to avoid the conversion to gtf/gff.
- Modifies `logger.py` to send logging information to stdout (this is also written to `postoga.log`)

# postoga 0.9.0-devel

- Adds `--outdir` to control where postoga output goes.
- Adds `--isoforms` to allow the user specify external isoform tables.
- Disables git features from `logger.py` [branch/commit] info.

# postoga 0.9.3-devel
- Re-implementation of postoga to match TOGA2.0 output
- Includes self-owned rust bed to gtf coverters through rustools
- Forces bed2gtf, bed2gff and gxf2bed installation for quick access
- Implements --engine option to use polars
- Now manages configuration and test through 'make configure' and 'make test'
- Adds license
- Drops --skip argument and adds --plot argument to plot stats [currently broken]
- Adds additional BUSCO and completeness stats to the main log file

# postoga 0.10.0-devel [BREAKING CHANGE]
- Replaces the ad-hoc logger with a shared `PostogaLogger`, adds a `-L/--level` flag (debug, info, warn, off), and pipes filter statistics directly into the logging stream.
- `--depure` now deletes every POSTOGA_* run artifact from the parent output directory before creating a new run, with warnings if files cannot be removed.
- Added `setup_postoga.sh` to bootstrap a `uv` virtual environment, install Python dependencies, and execute `maturin develop --release` for `rustools`.
- `filter_query_annotation()` now reports per-filter discard counts plus final projection/orthology summaries so CLI and `postoga.log` users can audit thresholds easily.
