// Bpipe pipeline to detect pathogenic STR expansions from exome data
// Set up to run on Meerkat cluster at MCRI

// Load system configuration and other settings
load 'pipeline_config.groovy'

// Load Bpipe pipeline stages
load 'pipeline_stages.groovy'

// Note: this will be populated by the set_sample_info stage
samples = Collections.synchronizedList([])

/////////////////////////////
// Run pipeline

run {
    '%_R*.fastq.gz' * [
        set_sample_info +
        align_bwa + index_bam +
        median_cov_region +
        STR_coverage +
        STR_locus_counts
    ] +
    estimate_size
}
