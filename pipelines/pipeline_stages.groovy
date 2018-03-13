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

    if(!file(REF).exists())
        fail """
             The configured decoy reference file: $REF could not be found.

             Please check pipelines/pipeline_config.groovy to make sure this is set correctly
        """

    [bwa,samtools,bedtools,goleft,python].each { tool ->
        if(!file(tool).exists())
            fail """
                 The location of tool $tool does not appear to exist.

                 Please check pipelines/pipeline_config.groovy to make sure this is set correctly
            """
    }

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
            set -o pipefail

            $bwa mem -M -t $threads
            -R "@RG\\tID:${sample}\\tPL:$PLATFORM\\tPU:NA\\tLB:${lane}\\tSM:${sample}"
            $REF
            $inputs |
            $samtools view -bSuh - | $samtools sort -o $output.bam -T $output.bam.prefix
        """, "bwamem"
    }
}

@preserve("*.bam")
align_bwa_bam = {
    doc "Align reads with bwa mem algorithm."

    def fastaname = get_fname(REF)
    produce(branch.sample + '.bam') {
        exec """
            set -o pipefail

            export JAVA_OPTS="-Dsamjdk.reference_fasta=$REF"

            gngstool ExtractFASTQ -bam ${input[input_type]} |
                $bwa mem -p -M -t $threads
                    -R "@RG\\tID:${sample}\\tPL:$PLATFORM\\tPU:NA\\tLB:${lane}\\tSM:${sample}"
                    $REF - |
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
            STRPATH=$PATH;
            PATH=$STRETCH/tools/bin:$PATH;
            $python $STRETCH/scripts/identify_locus.py
            --bam $input.bam
            --bed $STR_BED
            --output $output.locus_counts
            ;PATH=$STRPATH
        """
    }
}

estimate_size = {
    produce("STRs.tsv") {
        if(CONTROL=="") {
             exec """
                PATH=$PATH:$STRETCH/tools/bin;
                $python $STRETCH/scripts/estimateSTR.py
                    --locus_counts $inputs.locus_counts 
                    --STR_counts $inputs.STR_counts 
                    --median_cov $inputs.median_cov
                    --model $STRETCH/scripts/STRcov.model.csv
            """
        } else {
            exec """
                PATH=$PATH:$STRETCH/tools/bin;
                $python $STRETCH/scripts/estimateSTR.py
                    --locus_counts $inputs.locus_counts 
                    --STR_counts $inputs.STR_counts 
                    --median_cov $inputs.median_cov
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
        set -o pipefail

        $goleft covmed $input.bam | cut -f 1 > $output.median_cov
     """
}

/////////////////////////////////////
// Stages specific to exome pipeline

@transform('median_cov')
median_cov_region = {

doc "Calculate the median coverage over the target region"

    exec """
        set -o pipefail

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
            set -o pipefail

            $bedtools slop -b $SLOP -i $input.bed -g ${REF}.genome | $bedtools merge > $output.bed
        """
    //}
}

extract_reads_region = {

    doc "Extract reads from bam region + unaligned"

    def fastaname = get_fname(REF)

    produce(branch.sample + '_L001_R1.fastq.gz', branch.sample + '_L001_R2.fastq.gz') {
        exec """
            set -o pipefail

            cat <( $samtools view -hu -L $input.bed $input.bam )
                <( $samtools view -u -f 4 $input.bam ) |
            $samtools collate -Ou -n 128 - $output.prefix |
            $bedtools bamtofastq -i - -fq >(gzip -c > $output1.gz) -fq2 >(gzip -c > $output2.gz)
        """, "bwamem"
    }
}

@transform('median_cov')
median_cov_target = {

doc "Calculate the median coverage over the target region"

    exec """
        set -o pipefail

        $goleft covmed $input.bam $input.bed | cut -f 1 > $output.median_cov
     """
}
