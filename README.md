[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![version](https://img.shields.io/badge/version-0.4.0--devel-orange)

# postoga

The post-TOGA processing pipeline.


## What's new on version 0.4.0-devel

- splited main functions in ordered modules following TogaDir class layout (`filter_bed` will be improved in v.0.5.0-devel)
- improved logging (specifies some stats about annotation or filtered annotation [if user choose to filter it]) + connect() method to call it within a module (`bed2gtf`/`bed2gff` log info bothers the structure of postoga log -> will be fixed in v.0.5.0-devel)
- included `transcripts.quality.tsv` in query_table() and calculates some stats about it (included in logging)
- re-implemented `TOGA_assemblyStats.py` to estimate annotation completeness in single mode and haplotype-resolved assemblies -> `--assembly_qual` (to specify genes db -> tab-separated file [.txt,.tsv,.csv,...]; default ancestral); `--haplotype_path` (to specify path of both haplotypes) -> the equivalent of `TOGA_assemblyStats.py -m stats` && `TOGA_assemblyStats.py -m stats -ances ./ances.txt` is completed, `TOGA_assemblyStats.py -m merge` will be available in v.0.5.0-devel.
- postoga now automates: 1) dependencies [v.0.5.0-devel], 2) filtering, 3) conversion, 4) assembly quality estimation (just 1/4 quality steps postoga will have)
- shell is now `utils.py` (postoga base utility module)
- added https://github.com/hillerlab/TOGA/blob/master/TOGAInput/human_hg38/Ancestral_placental.txt to `./supply/Ancestral_placental.txt` 





## TO DO's

- [x] Build logger/version/constants scripts

- [x] postoga STEP 1: automate conversion from bed to gtf/gff -> bed2gtf will implement a sorting leaving dependecy on gtfsort (something already implemented in bed2gff

- [ ] Handle possible warnings of data (low number of intact genes, etc) within log file (?)

- [ ] Check/install/test requeriments (install_dependencies.py) -> rust, cargo, binaries (bed2gtf, bed2gff), compleasm (?)

- [x] postoga STEP 2: assembly stats re-implementation -> https://github.com/hillerlab/TOGA/blob/master/supply/TOGA_assemblyStats.py

- [ ] Build Graph class automating plotting reports -> need to define canvas structure and plot types to a default report (relevant variables)

- [ ] postoga STEP 3: alignment stats -> rust binding + ortholog lengths (?)

- [ ] postoga STEP 4: implement [compleasm](https://github.com/huangnengCSU/compleasm) an efficient alternative to BUSCO -> could be compared with assembly_stats.py -> then compared with an existente database of mammalian/placental genomes (?)

- [ ] ask Bogdan/Michael/Scott feedback and ideas