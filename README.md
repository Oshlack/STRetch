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
- Then activate the environment with `source activate STR`.
- For some cluster environments you may need to point the pipeline to the
specific version of Python, which would look something like this
`[where you installed conda]/miniconda3/envs/STR/bin/python`

R 3 (tested on R version 3.3.1)
- install.packages(c('optparse','plyr','dplyr','tidyr','reshape2'))

Command line tools
- BWA (version with bwa mem)
- Samtools

These commandline tools are required, but will already be installed by 
environment.yml as part of the python dependencies
- BedTools
- goleft

### Reference data

Required reference files:
- Reference genome with STR decoy chromosomes (fasta)
- BWA indices of reference
- .genome file of reference
- STR decoy bed file
- STR positions in genome annotated bed file

All chromosomes names in these files must be in the same sort order.

STRetch requires a reference genome that includes STR decoy chromosomes.
You can generate this by adding concatenating STRdecoys.fasta to the end of
any reference genome in fasta format, and then indexing the result.
For example:
```
cat ucsc.hg19.fasta STRdecoys.fasta > ucsc.hg19.STRdecoys.fasta
bwa index ucsc.hg19.STRdecoys.fasta
```

STRetch also requires a bed file specifying the position of all STRs in the
reference genome, with two additional columns containing the repeat unit/motif
and the number of repeat units in the reference.

For example:

```chr1	10000	10468	TAACCC	77.2```

This file is produced by extracting the appropriate columns from Tandem Repeats
Finder output.

You can download a bundle of pre-indexed genomes with STR decoys and
corresponding bed files here:
- [TODO]

### Configuring STRetch

Create a text file `STRetch/pipelines/pipeline_config.groovy` (note this file
must be in the same directory as the pipelines).
This will contain system-specific settings, especially paths to the locations
of software and reference files.
There is a template, and several examples provided in the form
`pipeline_config_[cluster-name].groovy`.


## Running STRetch

Run the appropriate pipeline (exome or WGS) using bpipe:
```bpipe run STRetch/pipelines/STRetch_exome_pipeline_meerkat.groovy sample_L001_R1.fastq.gz sample_L001_R2.fastq.gz```

All reads are assumed to be paired end with filenames ending in R1 and R2.
Multiple samples can be run in the same command and will be processed in
parallel and their STR variation compared to find outliers. For example:

```bpipe run STRetch/pipelines/STRetch_exome_pipeline_meerkat.groovy sample1_L001_R1.fastq.gz sample1_L001_R2.fastq.gz sample2_L001_R1.fastq.gz sample2_L001_R2.fastq.gz sample3_L001_R1.fastq.gz sample3_L001_R2.fastq.gz ...```

Note: for exome samples, the exome target region must be set in the
pipeline_config and is assumed to be the same for all samples.

## License

The code is freely available under the
[MIT license](http://www.opensource.org/licenses/mit-license.html).
