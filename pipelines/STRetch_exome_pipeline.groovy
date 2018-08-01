// Bpipe pipeline to detect pathogenic STR expansions from exome data
// Set up to run on Meerkat cluster at MCRI

// Load system configuration and other settings
load 'pipeline_config.groovy'

// Load Bpipe pipeline stages
load 'pipeline_stages.groovy'

/////////////////////////////
// Run pipeline

run {

    ~'(.*?)_.*_R[12].fastq.gz' * [
        set_sample_info +
        ~'.*?_(.*)_R[12].fastq.gz' * [
            set_fastq_info +
            align_bwa + index_bam
        ] + merge_bams +
        median_cov_region +
        STR_coverage +
        STR_locus_counts
    ] +
    estimate_size
}
