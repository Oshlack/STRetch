// Bpipe pipeline config file for:
// Bpipe pipeline to detect pathogenic STR expansions from exome data
// Set up to run on Meerkat cluster at MCRI

//STRetch installation location
STRETCH='..'

// Decoy reference assumed to have matching .genome file in the same directory
REF_DIR='$STRETCH/reference-data'
REF=REF_DIR + '/hg19.STRdecoys.sorted.fasta'
STR_BED=REF_DIR + '/hg19.simpleRepeat_period1-6.bed'
DECOY_BED=REF_DIR + '/STRdecoys.sorted.bed'

// For exome pipeline only
EXOME_TARGET="/group/bioi1/harrietd/ref-data/hg19_RefSeq_coding.sorted.bed"

// Software
PYTHON='/group/bioi1/harrietd/src/miniconda3/envs/STR/bin/python'

// Adjust parameters
PLATFORM='illumina'

// Number of threads to use for BWA
threads=8
