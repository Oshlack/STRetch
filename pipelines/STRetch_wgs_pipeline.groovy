// Bpipe pipeline to detect pathogenic STR expansions from whole genome sequencing data data
// Set up to run on the Broad cluster

// Load system configuration and other settings
load 'pipeline_config.groovy'

// Load Bpipe pipeline stages
load 'pipeline_stages.groovy'

run {
    '%_R*.fastq.gz' * [
        set_sample_info +
        align_bwa + index_bam +
        median_cov +
        STR_locus_counts 
    ] +
    estimate_size
}
