// Bpipe pipeline to extract paired end reads from a bam file

@filter('sorted')
sort_bam = {
// Sort by read name
    exec "samtools sort -n $input.bam $output.prefix"
}

extract_fastq = {
    def samplename = branch.name
    produce(samplename + '_L001_R1.fastq', samplename + '_L001_R2.fastq') {
        exec "bedtools bamtofastq -i $input.bam -fq $output1.fastq -fq2 $output2.fastq"
    }
}

gzip = {
    transform('.fastq') to('.fastq.gz') {
        exec "gzip -c $input.fastq > $output" 
    }
}

run {
    "%.bam" * [
        sort_bam + 
        extract_fastq
    ] + 
    "%.fastq" * [ 
        gzip 
    ]
}
