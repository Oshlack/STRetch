// Bpipe pipeline to detect pathogenic STR expansions from exome data
// Set up to run on Meerkat cluster at MCRI

// Load system configuration and other settings
load 'pipeline_config.groovy'

// Load Bpipe pipeline stages
load 'pipeline_stages.groovy'

/////////////////////////////
// Run pipeline

run {
    '%_R*.fastq.gz' * [
        set_sample_info +
        align_bwa + index_bam +
        median_cov +
        STR_coverage +
        STR_locus_counts +
        estimate_size
    ]
}
