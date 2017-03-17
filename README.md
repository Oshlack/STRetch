Note: this software is still very much in development, so use at own risk, and feel free to report any bugs.

# STRetch

Method for detecting pathogenic STR expansions from next-gen sequencing data

## Installing STRetch

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
- BWA (requires mem algorithm - tested on version 0.7.12)
- Samtools (tested with version 1.3.1, some older versions have different options so may require you to adapt the pipeline)

These commandline tools are required, but will already be installed by
environment.yml as part of the python dependencies
- BedTools
- goleft

### Reference data

**Quick start for human data:**
You can download a bundle of pre-indexed genomes with STR decoys and
corresponding bed files [here](https://figshare.com/s/1a39be9282c90c4860cd).

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
$ cat ucsc.hg19.fasta STRdecoys.fasta > ucsc.hg19.STRdecoys.fasta
$ bwa index ucsc.hg19.STRdecoys.fasta
```

STRetch also requires a bed file specifying the position of all STRs in the
reference genome, with two additional columns containing the repeat unit/motif
and the number of repeat units in the reference.

For example:

```chr1	10000	10468	TAACCC	77.2```

This file is produced by extracting the appropriate columns from Tandem Repeats
Finder output.

### Sample Installation

This is an example of installing STRetch in my home directory on the Broad
cluster, where BWA and Samtools, Bpipe, Conda and R are already installed:

```
$ cd /home/unix/hdashnow/git
$ git clone git@github.com:hdashnow/STRetch.git

# Download reference data (approx. 12 Gig, so could take a while)
$ cd /group/bioi1/harrietd/git/STRetch
$ wget -O reference-data.zip https://ndownloader.figshare.com/articles/4658701?private_link=1a39be9282c90c4860cd
$ unzip reference-data.zip

# Install conda environment
$ cd /home/unix/hdashnow/git/STRetch
$ conda env update -f environment.yml
# Install R packages
$ R
> install.packages(c('optparse','plyr','dplyr','tidyr','reshape2'))
> q()

# Check locations of software
$ source activate STR # Activates the Conda environment
$ which python
$ which goleft

# Set locations of software and reference files in pipeline_config.groovy
$ cd /home/unix/hdashnow/git/STRetch/pipelines
$ vim pipeline_config.groovy

//STRetch installation location
STRETCH='/home/unix/hdashnow/git/STRetch'
// Software
PYTHON='~/.conda/envs/STR/bin/python'
GOLEFT='~/.conda/envs/STR/bin/goleft'

$ cd /home/unix/hdashnow/git/STRetch
$ ln -s /home/unix/hdashnow/storage/ref-data/decoySTR ./reference-data
```

### Configuring STRetch

Create a text file `STRetch/pipelines/pipeline_config.groovy` (note this file
must be in the same directory as the pipelines).
This will contain system-specific settings, especially paths to the locations
of software and reference files.
There is a template, and several examples provided in the form
`pipeline_config_[cluster-name].groovy`.

### Sample Installation

This is an example of installing STRetch in my home directory on the Broad
cluster, where BWA and Samtools, Bpipe, Conda and R are already installed:

```
$ cd /home/unix/hdashnow/git
$ git clone git@github.com:hdashnow/STRetch.git
```
Download reference data (approx. 12 Gig, so could take a while, I find
downloading directly from Figshare using the browser to be faster)
```
$ cd /group/bioi1/harrietd/git/STRetch
$ mkdir reference-data
$ cd reference-data
$ wget -O reference-data.zip https://ndownloader.figshare.com/articles/4658701?private_link=1a39be9282c90c4860cd
$ unzip reference-data.zip
```
Install conda environment
```
$ cd /home/unix/hdashnow/git/STRetch
$ conda env update -f environment.yml
```
Install R packages
```
$ R
> install.packages(c('optparse','plyr','dplyr','tidyr','reshape2'))
> q()

```
Check locations of software
```
$ source activate STR # Activates the Conda environment
$ which python
$ which goleft
```
Set locations of software and reference files in pipeline_config.groovy
```
$ cd /home/unix/hdashnow/git/STRetch/pipelines
$ cp pipeline_config_template.groovy pipeline_config.groovy
```
Edit pipeline_config.groovy as appropriate, for example I needed to set:
```
//STRetch installation location
STRETCH='/home/unix/hdashnow/git/STRetch'
// Software
PYTHON='~/.conda/envs/STR/bin/python'
GOLEFT='~/.conda/envs/STR/bin/goleft'
```


## Running STRetch

Run the appropriate pipeline (exome or WGS) using bpipe:

```bpipe run STRetch/pipelines/STRetch_exome_pipeline_meerkat.groovy sample_L001_R1.fastq.gz sample_L001_R2.fastq.gz```

All reads are assumed to be paired end with filenames ending in
`_R1.fastq.gz` and `_R2.fastq.gz`.
Multiple samples can be run in the same command and will be processed in
parallel and their STR variation compared to find outliers. For example:

```bpipe run STRetch/pipelines/STRetch_exome_pipeline_meerkat.groovy sample1_L001_R1.fastq.gz sample1_L001_R2.fastq.gz sample2_L001_R1.fastq.gz sample2_L001_R2.fastq.gz sample3_L001_R1.fastq.gz sample3_L001_R2.fastq.gz ...```

Note: for exome samples, the exome target region must be set in the
pipeline_config and is assumed to be the same for all samples.

### Sample data

To quickly try out STRetch, you can
[download some simulated reads](https://figshare.com/s/cc7347f4637d9a7fe22d)
from samples with STR expansions in the SCA8 STR locus.

You will need to edit pipeline_config.groovy to use the provided target bed file:
```
EXOME_TARGET="SCA8_region.bed"
```

Sample command to run:
```
$ bpipe run path/to/STRetch/pipelines/STRetch_exome_pipeline.groovy *.fastq.gz
```

## License

The code is freely available under the
[MIT license](http://www.opensource.org/licenses/mit-license.html).
