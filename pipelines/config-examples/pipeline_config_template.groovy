// Bpipe pipeline config file
// Paths are relative to the directory the pipeline is running in, so absolute
// paths are recommended.

// Adjust parameters
PLATFORM='illumina'

// Number of threads to use for BWA
threads=8

// For exome pipeline only ***Edit before running the exome pipeline***
EXOME_TARGET=""

// STRetch installation location
STRETCH=""

// Paths to tools used by the pipeline
bpipe="$STRetch/tools/bin/bpipe"
python="$STRetch/tools/bin/python"
goleft="$STRetch/tools/bin/goleft"
bedtools="$STRetch/tools/bin/bedtools"
bwa="$STRetch/tools/bin/bwa"
samtools="$STRetch/tools/bin/samtools"

// Path to reference data
refdir="$STRetch/reference-data"

// Decoy reference assumed to have matching .genome file in the same directory
REF="$refdir/hg19.STRdecoys.sorted.fasta"
STR_BED="$refdir/hg19.simpleRepeat_period1-6_dedup.sorted.bed"
DECOY_BED="$refdir/STRdecoys.sorted.bed"

R="/usr/local/installed/R/3.3.2/lib64/R/bin/R"
