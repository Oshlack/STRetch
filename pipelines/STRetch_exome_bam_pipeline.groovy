// Bpipe pipeline to detect pathogenic STR expansions from exome sequencing data
// Takes a mapped bam as input and extracts relevant reads for remapping to STR decoys
// Usage:
// STRetch/tools/bin/bpipe run -p input_regions=all_STRs_in_genome.bed -p EXOME_TARGET=target_region.bed
// STRetch/pipelines/STRetch_wgs_bam_pipeline.groovy sample1.bam sample2.bam
// Note, EXOME_TARGET can also be set in pipeline_config.groovy

// Load system configuration and other settings
load 'pipeline_config.groovy'

// Load Bpipe pipeline stages
load 'pipeline_stages.groovy'

input_type='bam'

inputs "$input_type" : "Please supply one or more $input_type files to process"

bwa_parallelism = 1

shards = 1..bwa_parallelism

init_shard = {
    branch.shard = branch.name
}

run {
    "%.${input_type}" * [
        set_sample_info +
        mosdepth_region + mosdepth_mean_region +
        shards * [
            init_shard + align_bwa_bam + index_bam
        ] + merge_bams +
        STR_coverage +
        STR_locus_counts
    ] +
    estimate_size
}
