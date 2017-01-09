// Bpipe pipeline to detect pathogenic STR expansions from exome data

// Variables are set in pipeline_config.groovy

///////////////////
// Helper functions

def get_fname(path) {
    x = path.split('/')[-1]
    return(x)
}

//////////////////////////////////
// Stages common to all pipelines


/////////////////////////////////////
// Stages specific to exome pipeline


set_sample_info = {

    doc "Set information about the sample to be processed"

    branch.sample = branch.name

    }

@preserve("*.bam")
align_bwa = {
    doc "Concatenate with background reads then align with bwa mem algorithm."

    def fastaname = get_fname(REF)
    from('fastq.gz', 'fastq.gz') produce(branch.name + '.bam') {
        exec """
            bwa mem -M -t $threads
            -R "@RG\\tID:${sample}\\tPL:$PLATFORM\\tPU:1\\tLB:${sample}\\tSM:${sample}"
            $DECOY_REF
            $input1.gz
            $input2.gz |
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
            -g ${DECOY_REF}.genome
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

@transform('median_cov')
median_cov = {

doc "Calculate the median coverage over the target region"

    exec """
        goleft covmed $input.bam $EXOME_TARGET | cut -f 1 > $output.median_cov 
     """
}

estimate_size = {
        exec "$STRETCH/scripts/estimateSTR.R"
}

///////////////////////////////////
// Stages specific to WGS pipeline


