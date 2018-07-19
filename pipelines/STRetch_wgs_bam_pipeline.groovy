// Bpipe pipeline to detect pathogenic STR expansions from whole genome sequencing data
// Takes a mapped bam as input and extracts relevant reads for remapping to STR decoys

// Load system configuration and other settings
load 'pipeline_config.groovy'

// Load Bpipe pipeline stages
load 'pipeline_stages.groovy'

if(args.any { it.endsWith('.cram') })
    input_type = 'cram'
else
    input_type='bam'

inputs "$input_type" : "Please supply one or more $input_type files to process",

bwa_parallelism = 1

shards = 1..bwa_parallelism

if(input_type == "cram") 
    requires CRAM_REF: "To use CRAM format, please set the CRAM_REF parameter in pipeline_config.groovy to specify the reference to used to compress the CRAM file"

init_shard = {
    branch.shard = branch.name
}

run {
    "%.${input_type}" * [
        set_sample_info + 
        [ 
            mosdepth_dist + mosdepth_median,
            shards * [
                init_shard + align_bwa_bam + index_bam 
            ] + merge_bams
        ] +
        STR_coverage +
        STR_locus_counts 
    ] +
    estimate_size
}
