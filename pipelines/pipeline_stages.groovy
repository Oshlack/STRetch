// Bpipe pipeline to detect pathogenic STR expansions from exome data

// Variables are set in pipeline_config.groovy

///////////////////
// Helper functions

def get_fname(path) {
    x = path.split('/')[-1]
    return(x)
}

def get_info(filename) {
    return(filename.split('/')[-1].split('\\.')[0].split('_'))
}

//////////////////////////////////
// Stages common to all pipelines

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    def info = get_info(input)
    branch.sample = info[0]
    branch.lane = info[1]
}

@preserve("*.bam")
align_bwa = {
    doc "Align reads with bwa mem algorithm."

    def fastaname = get_fname(REF)
    from('fastq.gz', 'fastq.gz') produce(branch.sample + '.bam') {
        exec """
            bwa mem -M -t $threads
            -R "@RG\\tID:${sample}\\tPL:$PLATFORM\\tPU:NA\\tLB:${lane}\\tSM:${sample}"
            $REF
            $inputs |
            samtools view -bSuh - | samtools sort -o $output.bam -T $output.bam.prefix
        """, "bwamem"
    }
}

@preserve("*.bai")
index_bam = {
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam"
    }
    forward input
}

STR_coverage = {
    transform("bam") to ("STR_counts") {
        exec """
            bedtools coverage -counts
            -sorted
            -g ${REF}.genome
            -a $DECOY_BED
            -b $input.bam > $output.STR_counts
        """
    }
}

STR_locus_counts = {
    transform("bam") to ("locus_counts") {
        exec """
            $PYTHON $STRETCH/scripts/identify_locus.py
            --bam $input.bam
            --bed $STR_BED
            --output $output.locus_counts
        """
    }
}

estimate_size = {
        exec "$STRETCH/scripts/estimateSTR.R"
}

///////////////////////////////////
// Stages specific to WGS pipeline

@transform('median_cov')
median_cov = {

doc "Calculate the median coverage over the whole genome"

    exec """
        goleft covmed $input.bam | cut -f 1 > $output.median_cov
     """
}

/////////////////////////////////////
// Stages specific to exome pipeline

@transform('median_cov')
median_cov_region = {

doc "Calculate the median coverage over the target region"

    exec """
        goleft covmed $input.bam $EXOME_TARGET | cut -f 1 > $output.median_cov
     """
}
