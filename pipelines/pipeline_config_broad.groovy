// Bpipe pipeline config file for:
// Bpipe pipeline to detect pathogenic STR expansions from whole genome sequencing data data
// Set up to run on the Broad cluster

//STRetch installation location
STRETCH='..'

REF_DIR='$STRETCH/reference-data'
REF=REF_DIR + '/hg19.STRdecoys.fasta'
DECOY=REF_DIR + '/STRdecoys.fasta'
DECOY_BED=REF_DIR + '/STRdecoys.bed'
STR_BED=REF_DIR + '/hg19.simpleRepeat_period1-6.bed'
PYTHON='~/.conda/envs/STR/bin/python'
PICARD='/seq/software/picard/current/bin/picard.jar'

// Adjust parameters
PLATFORM='illumina'

// Number of threads to use for BWA
threads=8
