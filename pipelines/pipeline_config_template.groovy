// Template Bpipe pipeline config file
// Paths are relative to the directory the pipeline is running in, so absolute
// paths are recommended.

//STRetch installation location
STRETCH='/absolute/path/to/installation/STRetch'

// Decoy reference assumed to have matching .genome file in the same directory
DECOY_REF="$STRETCH/reference-data/hg19.STRdecoys.sorted.fasta"
STR_BED="$STRETCH/reference-data/hg19.simpleRepeat_period1-6.bed"
DECOY_BED="$STRETCH/reference-data/STRdecoys.sorted.bed"

// For exome pipeline only
EXOME_TARGET="path/to/exome_target_regions.bed"

// Software
// you may need to set absolute paths for these if source active STR is no feasible
// e.g. python would be found in '/path/to/miniconda3/envs/STR/bin/python'
PYTHON='python'
GOLEFT='goleft'

// Adjust parameters
PLATFORM='illumina'

// Number of threads to use for BWA
threads=8
