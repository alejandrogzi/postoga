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
- postoga now automates installing requirements (python/rust)