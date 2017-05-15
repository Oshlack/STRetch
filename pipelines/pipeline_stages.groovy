// Bpipe pipeline to detect pathogenic STR expansions from exome data

// Variables e.g. tools are set in pipeline_config.groovy

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
    if (info.length >= 2) {
        branch.lane = info[1]
    } else {
        branch.lane = 'L001'
    }

}

@preserve("*.bam")
align_bwa = {
    doc "Align reads with bwa mem algorithm."

    def fastaname = get_fname(REF)
    from('fastq.gz', 'fastq.gz') produce(branch.sample + '.bam') {
        exec """
            $bwa mem -M -t $threads
            -R "@RG\\tID:${sample}\\tPL:$PLATFORM\\tPU:NA\\tLB:${lane}\\tSM:${sample}"
            $REF
            $inputs |
            $samtools view -bSuh - | $samtools sort -o $output.bam -T $output.bam.prefix
        """, "bwamem"
    }
}

@preserve("*.bai")
index_bam = {
    transform("bam") to ("bam.bai") {
        exec "$samtools index $input.bam"
    }
    forward input
}

STR_coverage = {
    transform("bam") to ("STR_counts") {
        exec """
            $bedtools coverage -counts
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
            source $STRETCH/tools/bin/activate STR; python $STRETCH/scripts/identify_locus.py
            --bam $input.bam
            --bed $STR_BED
            --output $output.locus_counts
        """
    }
}

estimate_size = {
    produce("STRs.tsv") {
        if(CONTROL=="") {
             exec """
                Rscript $STRETCH/scripts/estimateSTR.R 
                    --model $STRETCH/scripts/STRcov.model.csv 
            """
        } else {
            exec """
                Rscript $STRETCH/scripts/estimateSTR.R 
                    --model $STRETCH/scripts/STRcov.model.csv 
                    --control $CONTROL
            """
        }
    }
}

///////////////////////////////////
// Stages specific to WGS pipeline

@transform('median_cov')
median_cov = {

doc "Calculate the median coverage over the whole genome"

    exec """
        $goleft covmed $input.bam | cut -f 1 > $output.median_cov
     """
}

/////////////////////////////////////
// Stages specific to exome pipeline

@transform('median_cov')
median_cov_region = {

doc "Calculate the median coverage over the target region"

    exec """
        $goleft covmed $input.bam $EXOME_TARGET | cut -f 1 > $output.median_cov
     """
}

///////////////////////////////////
// Stages specific to WGS bam pipeline

@filter('slop')
str_targets = {
        
    doc "Create bed file of region likely to contain STR reads and their mates"

    SLOP=800

    //produce(STR_BED[0..-3] + 'slop.bed') {    
        exec """
            $bedtools slop -b $SLOP -i $input.bed -g ${REF}.genome | $bedtools merge > $output.bed
        """
    //}
}

extract_reads_region = {

    doc "Extract reads from bam region + unaligned"

    def fastaname = get_fname(REF)

    produce(branch.sample + '_L001_R1.fastq.gz', branch.sample + '_L001_R2.fastq.gz') {
        exec """

            cat <( $samtools view -hu -L $input.bed $input.bam ) 
                <( $samtools view -u -f 4 $input.bam ) | 
            $samtools collate -Ou -n 128 - $output.prefix | 
            $bedtools bamtofastq -i - -fq >(gzip -c > $output1.gz) -fq2 >(gzip -c > $output2.gz)
        """
    }
}

@transform('median_cov')
median_cov_target = {

doc "Calculate the median coverage over the target region"

    exec """
        $goleft covmed $input.bam $input.bed | cut -f 1 > $output.median_cov
     """
}
