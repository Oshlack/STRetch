// Bpipe pipeline to detect pathogenic STR expansions from exome data

// Variables e.g. tools are set in pipeline_config.groovy

// Amount to inflate STR regions by for remapping
SLOP=5

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

    branch.original_bam = input[input_type]

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

    var input_regions : false,
        pairThreads : 2

    String regionFlag = ""
    if(input_regions)
        regionFlag="-L $input_regions"


    String shardFlag = bwa_parallelism > 1 ? "-s $shard,$bwa_parallelism" : ""

    String fastaname = get_fname(REF)


    // :$GNGS_JAR 

    // we construct the output file name by joining the following parts:
    // - the sample id
    // - the shard number (only if parallelism in use)
    // - the file extension '.bam'
    List outputFileParts = [branch.sample] + 
                           (bwa_parallelism>1?[shard]:[]) + // empty unless parallelism used
                           ['bam']

    produce(outputFileParts.join('.')) {
        exec """
            set -o pipefail

            export JAVA_OPTS="-Dsamjdk.reference_fasta=$CRAM_REF"

            java -Xmx16g -Dsamjdk.reference_fasta=$CRAM_REF 
                 -jar $STRETCH/tools/bazam/build/libs/bazam.jar
                 -pad $SLOP -n 6
                $regionFlag $shardFlag -bam ${input[input_type]} |
                $bwa mem -p -M -t ${threads-1}
                    -R "@RG\\tID:${sample}\\tPL:$PLATFORM\\tPU:NA\\tLB:${lane}\\tSM:${sample}"
                    $REF - |
                $samtools view -bSuh - | $samtools sort -o $output.bam -T $output.bam.prefix
        """, "bwamem"
    }
}

merge_bams = {

    if(bwa_parallelism == 1) {
        println "Skipping merge because bwa parallelism not in use"
        return
    }

    produce(branch.sample + '.merge.bam') {
        exec """
            time java -Xmx2g -jar $STRETCH/tools/picard.jar MergeSamFiles
                ${inputs.bam.withFlag("INPUT=")}
                VALIDATION_STRINGENCY=LENIENT
                ASSUME_SORTED=true
                CREATE_INDEX=true
                OUTPUT=$output.bam
         """, "merge"

/*
        exec """
            $samtools merge $output.bam $inputs.bam
        """
*/

    }

}


@preserve("*.bai")
index_bam = {
    transform("bam") to("bam.bai") {
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

    from(original_bam) {
        exec """
            set -o pipefail

            $goleft covmed $input.bam | cut -f 1 > $output.median_cov
         """
    }
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


    //produce(STR_BED[0..-3] + 'slop.bed') {
        exec """
            set -o pipefail

            $bedtools merge -i $input.bed > $output.bed
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

mosdepth_dist = {
    
    doc "Calculate cumulative depth distribution from a bam or cram file"

    transform(input_type) to(input.prefix + '.mosdepth.global.dist.txt') {

        if(input_type=="cram") {
            exec """
                $STRETCH/tools/bin/mosdepth -n -t $threads
                -f $CRAM_REF
                $output.prefix.prefix.prefix.prefix
                ${input[input_type]}
            ""","mosdepth"
        } else {
            exec """
                $STRETCH/tools/bin/mosdepth -n -t $threads
                $output.prefix.prefix.prefix.prefix
                ${input[input_type]}
            ""","mosdepth"
        }

    }
}

mosdepth_median = {
    transform('mosdepth.global.dist.txt') to('median_cov') {
        doc "Calculate the median coverage from mosdepth .dist.txt output"
    
        from('.dist.txt') {
            exec """
                $python $STRETCH/scripts/mosdepth_median.py --out $output.median_cov $input.txt
            """
        }
    }
}
