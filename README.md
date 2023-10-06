[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![version](https://img.shields.io/badge/version-0.2.0--devel-orange)

# postoga

The post-TOGA processing pipeline.


TO DO's

[x] Build logger/version/constants scripts

[x] postoga STEP 1: automate conversion from bed to gtf/gff -> bed2gtf will implement a sorting leaving dependecy on gtfsort (something already implemented in bed2gff)

[ ] Check/install/test requeriments (install_dependencies.py) -> rust, cargo, binaries (bed2gtf, bed2gff), compleasm (?)

[ ] postoga STEP 2: assembly stats re-implementation -> https://github.com/hillerlab/TOGA/blob/master/supply/TOGA_assemblyStats.py

[ ] Build Graph class automating plotting reports -> need to define canvas structure and plot types to a default report (relevant variables)

[ ] postoga STEP 3: alignment stats -> rust binding + ortholog lengths (?)

[ ] postoga STEP 4: implement [compleasm](https://github.com/huangnengCSU/compleasm) an efficient alternative to BUSCO -> could be compared with assembly_stats.py -> then compared with an existente database of mammalian/placental genomes (?)

[ ] ask Bogdan/Michael/Scott feedback and ideas