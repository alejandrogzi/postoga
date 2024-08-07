> [!NOTE]
>
> postoga is dependent from [TOGA](https://github.com/hillerlab/TOGA). Any changes in TOGA will have a repercusion here. If you found any bug/errors, please report them here.
> This project is in constant development, any desired features are welcome!


# postoga

The post-[TOGA](https://github.com/hillerlab/TOGA) processing pipeline.

![version](https://img.shields.io/badge/version-0.9.2--devel-orange)

<p align="center">
    <img width=700 align="center" src="./supply/postoga_logo_git.png" >
</p>

<!-- <img src="./supply/postoga_report.png" align="center"/> -->

> ## What's new on version 0.9.2-devel
>
> - CI workflow + testing
> - Enhanced env creation (using hatch)

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
pip install hatch && hatch shell
hatch run configure
```

Run test to confirm functionality
```bash
hatch run test
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

> [!TIP]
>
> If the only thing you want to do is apply some filters to a TOGA result or convert results to GTF/GFF files, I recommend the following command:
>
> ```
> ./postoga.py base \
> --togadir /your/TOGA/dir \
> --outdir /your/out/dir \ 
> -bc [YOUR CLASSES] \
> -br [YOUR RELATIONS] \
> -th [YOUR THRESHOLD] \
> -to [YOUR FORMAT GTF/GFF/BED] \
> --skip
> ```

```text
usage: postoga.py [-h] {base,haplotype}

positional arguments:
  {base,haplotype}  Select mode
    base            Base mode
    haplotype       Haplotype mode

postoga.py base [-h] --outdir OUTDIR --togadir TOGADIR [-bc BY_CLASS] [-br BY_REL]
                       [-th THRESHOLD] -to {gtf,gff,bed} [-aq ASSEMBLY_QUAL]
                       [-sp {human,mouse,chicken}] [-src {ensembl,gene_name,entrez}]
                       [-phy {mammals,birds}] [-s] [-par PARALOG]

optional arguments:
  -h, --help            show this help message and exit
  --outdir OUTDIR, -o OUTDIR
                        Path to posTOGA output directory
  --togadir TOGADIR, --td TOGADIR
                        Path to TOGA results directory
  -bc BY_CLASS, --by-class BY_CLASS
                        Filter parameter to only include certain orthology classes (I, PI, UL, M, PM,
                        L, UL)
  -br BY_REL, --by-rel BY_REL
                        Filter parameter to only include certain orthology relationships (o2o, o2m,
                        m2m, m2m, o2z)
  -th THRESHOLD, --threshold THRESHOLD
                        Filter parameter to preserve orthology scores greater or equal to a given
                        threshold (0.0 - 1.0)
  -to {gtf,gff,bed}, --to {gtf,gff,bed}
                        Specify the conversion format for .bed (query_annotation/filtered) file (gtf,
                        gff3) or just keep it as .bed (bed)
  -aq ASSEMBLY_QUAL, --assembly_qual ASSEMBLY_QUAL
                        Calculate assembly quality based on a list of genes provided by the user
                        (default: Ancestral_placental.txt)
  -sp {human,mouse,chicken}, --species {human,mouse,chicken}
                        Species name to be used as a reference for the assembly quality calculation
                        (default: human)
  -src {ensembl,gene_name,entrez}, --source {ensembl,gene_name,entrez}
                        Source of the ancestral gene names (default: ENSG)
  -phy {mammals,birds}, --phylo {mammals,birds}
                        Phylogenetic group of your species (default: mammals)
  -s, --skip            Skip steps 2, 3, and 4 and only filter the .bed file
  -par PARALOG, --paralog PARALOG
                        Filter parameter to preserve transcripts with paralog projection probabilities
                        less or equal to a given threshold (0.0 - 1.0)
  -iso ISOFORMS, --isoforms ISOFORMS
                        Path to a custom isoform table (default: None)

postoga.py haplotype [-h] --outdir OUTDIR -hp HAPLOTYPE_DIR [-r RULE] [-s {query,loss}]

optional arguments:
  -h, --help            show this help message and exit
  --outdir OUTDIR, -o OUTDIR
                        Path to posTOGA output directory
  -hp HAPLOTYPE_DIR, --haplotype_dir HAPLOTYPE_DIR
                        Path to TOGA results directories separated by commas (path1,path2,path3)
  -r RULE, --rule RULE  Rule to merge haplotype assemblies (default: I>PI>UL>L>M>PM>PG>abs)
  -s {query,loss}, --source {query,loss}
                        Source of the haplotype classes (query, loss)
```


