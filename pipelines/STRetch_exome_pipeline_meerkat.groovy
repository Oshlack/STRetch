// Bpipe pipeline to detect pathogenic STR expansions from exome data
// Set up to run on Meerkat cluster at MCRI

//STRetch installation location
STRETCH='../'

// Decoy reference assumed to have matching .genome file in the same directory
DECOY_REF="$STRETCH/reference-data/hg19.STRdecoys.sorted.fasta"
EXOME_TARGET="$STRETCH/reference-data/hg19_RefSeq_coding.sorted.bed"
STR_BED=REF_DIR + '/hg19.simpleRepeat_period1-6.bed'

// Software
PYTHON='/group/bioi1/harrietd/src/miniconda3/envs/STR/bin/python'

// Adjust parameters
PLATFORM='illumina'

def get_fname(path) {
    x = path.split('/')[-1]
    return(x)
}

/////////////////////////////
// Stages

set_sample_info = {

    doc "Set information about the sample to be processed"

    branch.sample = branch.name

    }

threads=8

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
            -a $STRETCH/reference-data/STRdecoys.bed
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

@transform('coverage')
coverage = {
    exec """
        bedtools coverage
            -sorted -d
            -g ${DECOY_REF}.genome
            -a $EXOME_TARGET
            -b $input.bam
            > $output.coverage
     ""","bedtools"
}

@transform('median_cov')
median_cov = {
    exec """
        cut -f 8 $input.coverage  | sort -n | awk -f $STRETCH/scripts/median.awk > $output.median_cov
     """
}

estimate_size = {
        exec "$STRETCH/scripts/estimateSTR.R"
}

/////////////////////////////
// Run pipeline


run {
    '%_R*.fastq.gz' * [
        set_sample_info +
        align_bwa + index_bam +
        STR_locus_counts +
        coverage +
        median_cov +
        STR_coverage +
        estimate_size
    ]
}
