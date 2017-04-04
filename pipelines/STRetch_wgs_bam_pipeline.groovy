// Bpipe pipeline to detect pathogenic STR expansions from whole genome sequencing data
// Takes a mapped bam as input and extracts relevant reads for remapping to STR decoys

// Load system configuration and other settings
load 'pipeline_config.groovy'

// Load Bpipe pipeline stages
load 'pipeline_stages.groovy'

run {
    str_targets +
    '%.bam' * [
        set_sample_info +
        extract_reads_region +
        align_bwa + index_bam +
        median_cov_target +
        STR_coverage +
        STR_locus_counts 
    ] +
    estimate_size
}
