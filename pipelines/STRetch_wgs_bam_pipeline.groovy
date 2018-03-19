// Bpipe pipeline to detect pathogenic STR expansions from whole genome sequencing data
// Takes a mapped bam as input and extracts relevant reads for remapping to STR decoys

// Load system configuration and other settings
load 'pipeline_config.groovy'

// Load Bpipe pipeline stages
load 'pipeline_stages.groovy'

input_type='bam'

inputs "$input_type" : "Please supply one or more $input_type files to process",
       "bed"         : "Please give a BED file defining the target regions to analyse"

bwa_parallelism = 1

shards = 1..bwa_parallelism

run {
    str_targets +
    "%.${input_type}" * [
        set_sample_info +
        shards * [
            { branch.shard = branch.name } + align_bwa_bam + index_bam 
        ] + merge_bams +
        median_cov_target +
        STR_coverage +
        STR_locus_counts 
    ] +
    estimate_size
}
