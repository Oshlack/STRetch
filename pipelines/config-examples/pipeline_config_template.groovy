// Bpipe pipeline config file
// Paths are relative to the directory the pipeline is running in, so absolute
// paths are recommended.

// Adjust parameters
PLATFORM='illumina'

// Number of threads to use for BWA
threads=8

// For exome pipeline only ***Edit before running the exome pipeline***
EXOME_TARGET="path/to/exome_target_regions.bed" 
