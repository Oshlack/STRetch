// Bpipe pipeline to detect pathogenic STR expansions from exome data
// Set up to run on the Broad cluster

//STRetch installation location
STRETCH='../'

REF_DIR='$STRETCH/reference-data'
REF=REF_DIR + '/hg19.STRdecoys.fasta'
DECOY=REF_DIR + '/STRdecoys.fasta'
DECOY_BED=REF_DIR + '/STRdecoys.bed'
STR_BED=REF_DIR + '/hg19.simpleRepeat_period1-6.bed'
PYTHON='~/.conda/envs/STR/bin/python'
PICARD='/seq/software/picard/current/bin/picard.jar'

def get_info(filename) {
    return(filename.split('/')[-1].split('\\.')[0].split('_'))
}

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    def info = get_info(input)
    branch.sample = info[0]
    branch.lane = info[1]
    }

PLATFORM='illumina'
threads=8

align_ref = {
    doc "Extract reads from bam then align with bwa mem algorithm."

    from('fastq.gz', 'fastq.gz') transform('.bam') {
//    from('fastq', 'fastq') transform('.bam') {
        exec """
            bwa mem -M
            -t $threads
            -R "@RG\\tID:${sample}_${lane}\\tPL:$PLATFORM\\tPU:NA\\tLB:${lane}\\tSM:${sample}"
            $REF $inputs |
            samtools view -bSh - > $output.bam
        """, "bwa"
    }
}

align_decoy = {
    doc "Extract reads from bam then align with bwa mem algorithm."

    from('fastq.gz', 'fastq.gz') transform('.decoy.bam') {
        exec """
            bwa mem -M
            -R "@RG\\tID:${sample}_${lane}\\tPL:$PLATFORM\\tPU:NA\\tLB:${lane}\\tSM:${sample}"
            $DECOY $inputs.gz |
            samtools view -bSuh - |
            samtools sort -o $output.bam -
        """, "bwa"
    }
}

@filter('sorted')
sort_bam = {
    exec "samtools sort $input.bam $output.bam.prefix", "medium"
}

index_bam = {
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam", "medium"
    }
    forward input
}

@transform('count')
count_bam = {
    exec """
        samtools view -c $input.bam > $output.count
     ""","medium"
}

median_cov_samtools = {
    exec """
        time samtools depth -a $input.bam -b $input.bed | cut -f 3  | sort -n | awk -f $STRETCH/scripts/median.awk > $output
     ""","medium"
}

@transform('coverage')
coverage = {
    exec """
        bedtools coverage
            -sorted -d
            -g ${REF}.genome
            -a $input.bed
            -b $input.bam
            > $output.coverage
     ""","medium"
}

median_cov = {
    exec """
        cut -f 5 $input.coverage  | sort -n | awk -f $STRETCH/scripts/median.awk > $output
     ""","medium"
}

STR_coverage_ref = {
    transform("bam") to ("STRcoverage.txt") {
        exec """
            bedtools coverage -counts
            -sorted
            -g ${REF}.genome
            -a $DECOY_BED
            -b $input.bam > $output.txt
        """
    }
}

STR_coverage_decoy = {
    transform("bam") to ("STRcoverage.txt") {
        exec """
            bedtools coverage -counts
            -sorted
            -g ${DECOY}.genome
            -a $DECOY_BED
            -b $input.bam > $output.txt
        """
    }
}

STR_locus_counts = {
    transform("bam") to ("locus_counts.txt") {
        exec """
            $PYTHON $STRETCH/scripts/identify_locus.py
            --bam $input.bam
            --bed $STR_BED
            --output $output.txt
        ""","medium"
    }
}

@transform('metrics')
metrics = {
    exec """
        java -Xmx60g -jar $PICARD CollectWgsMetrics
            R=$REF
            I=$input.bam
            O=$output.metrics
    ""","large"
}

insertsize = {
    exec """
        java -Xmx8g -jar $PICARD CollectInsertSizeMetrics
            HISTOGRAM_FILE=$output.pdf
            I=$input.bam
            O=$output.txt
    ""","medium"
}

estimate_size = {
        exec "$STRETCH/scripts/estimateSTR.R"
}

run {
    '%_R*.fastq.gz' * [
        set_sample_info +
        align_ref +
        sort_bam +
        index_bam +
        STR_coverage_ref +
        STR_locus_counts +
        metrics +
        insertsize +
        estimate_size
    ]
}
