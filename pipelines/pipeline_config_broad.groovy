// Bpipe pipeline config file for:
// Bpipe pipeline to detect pathogenic STR expansions from whole genome sequencing data data
// Set up to run on the Broad cluster

//STRetch installation location
STRETCH='..'

// Decoy reference assumed to have matching .genome file in the same directory
REF_DIR='$STRETCH/reference-data'
REF=REF_DIR + '/hg19.STRdecoys.fasta'
STR_BED=REF_DIR + '/hg19.simpleRepeat_period1-6.bed'
DECOY_BED=REF_DIR + '/STRdecoys.bed'

// Software
PYTHON='~/.conda/envs/STR/bin/python'

// Adjust parameters
PLATFORM='illumina'

// Number of threads to use for BWA
threads=8
