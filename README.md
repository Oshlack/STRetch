# STRetch

Method for detecting pathogenic STR expansions from next-gen sequencing data

## Installing STRetch
- [TODO]

### Requirements

[Bpipe](http://docs.bpipe.org/)

Python 3
- See `environment.yml` for specific packages
- Easiest to install using [Conda](http://conda.pydata.org/docs/using/envs.html):
`conda env create -f environment.yml`

R 3 (tested on R version 3.3.1)
- install.packages(c('optparse','plyr','dplyr','tidyr','reshape2'))

Command line tools
- BWA (version with bwa mem)
- Samtools
- BedTools
- Picard Tools

### Reference data

STRetch requires a reference genome that includes STR decoy chromosomes.
You can generate this by adding concatenating STRdecoys.fasta to the end of
any reference genome in fasta format, and then indexing the result.
For example:
```
cat ucsc.hg19.fasta STRdecoys.fasta > ucsc.hg19.STRdecoys.fasta
bwa index ucsc.hg19.STRdecoys.fasta
```

You can download pre-indexed genomes with STR decoys here:
- [TODO]

## Running STRetch
- [TODO]

## License

The code is freely available under the
[MIT license](http://www.opensource.org/licenses/mit-license.html).
