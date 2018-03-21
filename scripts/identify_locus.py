#!/usr/bin/env python
"""Look at reads mapping to decoy STR chromosomes and use their pairs to
determine which STR locus they originated from.
Generates counts of reads mapping to the decoy originating each STR locus.
"""

import argparse
import pysam
#from insert_size import parse_bed, span_intervals
import pybedtools as bt
import pandas as pd
import numpy as np
import os
import random
import string
from decoy_STR import normalise_str
import sys #XXX Just for testing

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Reports counts of reads mapping to the STR decoy originating each STR locus.')
    parser.add_argument(
        '--bam', type=str, required=True, nargs='+',
        help='One or more BAM files of mapped reads e.g. produced by Bowtie2 or BWA. If multiple BAMs all are assumed to come from the same sample.')
    parser.add_argument(
        '--bed', type=str, required=False,
        help='bed file containing genomic locations of STR loci. Genomic locations should be relative to the fasta reference used to create the bam')
    parser.add_argument(
        '--output', type=str, required=False,
        help='Output file name. Defaults to stdout.')
    parser.add_argument(
        '--dist', type=int, required=False, default=500,
        help='Counts are only assigned to an STR locus that is at most this many bp away. The best choice is probably the insert size. (default: %(default)s)')
    return parser.parse_args()

def randomletters(length):
   return ''.join(random.choice(string.ascii_lowercase) for i in range(length))

def detect_readlen(bamfile, sample = 100):
    """Samples reads from a bam file to detect the read length. Assumes uniform
    read lengths and that any shorter reads are due to trimming or clipping.
    Args:
        bamfile (str): Location of bamfile (will be opened and closed, so don't
            use an existing open bamfile handle).
        sample (int): Number of reads to sample.

    Returns:
        int: The maximum read length observed.
    """
    bam = pysam.Samfile(bamfile, 'rb')
    maxlen = 0
    count = 0
    count_noCIGAR = 0

    for read in bam.fetch():
        count += 1
        readlen = read.infer_query_length(always=False)
        if readlen:
            if readlen > maxlen:
                maxlen = readlen
        else:
            count_noCIGAR += 1

        if count >= sample:
            bam.close()
            return maxlen, count_noCIGAR

def main():
    # Parse command line arguments
    args = parse_args()
    bamfiles = args.bam
    bedfile = args.bed
    outfile = args.output

    max_distance = args.dist

    # Check bamfiles have unique names
    if len(set(bamfiles)) < len(bamfiles):
        sys.exit('ERROR: There were multiple bamfiles with the same filename. Please check your input')

    #STR_bed = parse_bed(args.bed, position_base=0)
    STR_bed = bt.BedTool(bedfile)

    all_results = []
    for bamfile in bamfiles:

        readlen, count_noCIGAR = detect_readlen(bamfile)

        # Read bam
        bam = pysam.Samfile(bamfile, 'rb')

        # Get relevant chromosomes
        required_chroms = []
        unpaired = 0
        total = 0
        for chrom in bam.references:
            if chrom.startswith('STR-'):
                required_chroms.append(chrom)
        # Check if any STR- chromosomes
        if len(required_chroms) == 0:
           sys.exit('ERROR: There were no reads mapping to chromosomes with names starting with "STR-" in {0}. Are you sure this data is mapped to a reference genome with STR decoy chromosomes?'.format(bamfile))

        for chrom in required_chroms:
            motif = chrom.split('-')[1]
            all_positions = []
            all_segments = bam.fetch(reference=chrom)

            for read in all_segments:
                total += 1
                try:
                    mate_chr = read.next_reference_name
                except ValueError:
                    unpaired += 1
                    continue
                mate_start = read.next_reference_start
                mate_stop = mate_start + readlen
                all_positions.append([mate_chr, mate_start, mate_stop])

            # Strategy:
            # Merge all overlapping intervals
            # Keep the count of reads corresponding to each merged interval (i.e. 1 for each read contained in it)
            # Assign each interval to the closest STR (within 500 bp? - the insert size) with the correct motif, adding together the count of reads
            # Check motif: == normalise_str()
            # There should be two 1-2 intervals per STR, likely one for each flank.
            # Report the read count for each STR

            if len(all_positions) > 0:
                motif_bed = bt.BedTool(all_positions).sort()
                # Merge all the intervals, then count how many of the original intervals overlap the merged ones (4th column)
                motif_coverage = motif_bed.merge(stream=True).coverage(b=motif_bed, counts=True)

                tmp_bed = 'tmp-' + randomletters(8) + '.bed' #create temporary file for bedtools to write to and pandas to read since streams don't seem to work
                closest_STR = motif_coverage.closest(STR_bed, d=True, stream=True).saveas(tmp_bed)
                colnames = ['chr', 'start', 'stop', 'count', 'STR_chr', 'STR_start',
                'STR_stop', 'motif', 'reflen', 'distance']
                df = pd.read_csv(tmp_bed, sep='\t', header=None, names=colnames)
                os.remove(tmp_bed) #delete temporary file

                # Filter out loci that are too far away
                df = df.loc[df['distance'] <= max_distance, :]
                df['motif'] = df['motif'].map(normalise_str) # Normalise the STR motif to enable comparisons
                # Remove STRs that don't match the motif
                df = df.loc[df['motif'] == normalise_str(motif), :]
                df = df.loc[:, ['STR_chr', 'STR_start', 'STR_stop', 'motif', 'count', 'reflen']]

                all_results.append(df)
        # Print a warning message in case of reads without a CIGAR string
        if count_noCIGAR > 0:
            sys.stderr.write('WARNING: ' + str(count_noCIGAR) + ' read(s) in ' + bamfile + ' file had no CIGAR string.')

        if total == 0:
            sys.exit('ERROR: there were no reads overlapping the target STR regions. This may indicate a problem with the input file.\n')
        elif unpaired == total:
            sys.exit('ERROR: all {0} reads overlapping the target STR regions appear to be unpaired. You may wish to check your bam file is paired-end and correctly formed.\n'.format(total))
        elif unpaired > 0:
            sys.stderr.write('WARNING: it appears that {0} of the {1} reads overlapping the target STR regions were unpaired and so no useful data could be obtained from them.\n'.format(unpaired, total))



    # Sum counts from multiple bam files and multiple rows
    if len(all_results) == 1:
        df_total = all_results[0]
    else:
        df_total = pd.concat(all_results, ignore_index=True)

    summed = df_total.groupby(['STR_chr', 'STR_start', 'STR_stop', 'motif', 'reflen'], as_index=False).aggregate(np.sum)

    # Print a warning message in case of reads without a CIGAR string
    if count_noCIGAR > 0:
        sys.stderr.write('WARNING: ' + str(count_noCIGAR) + ' read(s) in ' + bamfile + ' file had no CIGAR string.')

    if total == 0:
        sys.exit('ERROR: there were no reads overlapping the target STR regions. This may indicate a problem with the input file.\n')
    elif unpaired == total:
        sys.exit('ERROR: all {0} reads overlapping the target STR regions appear to be unpaired. You may wish to check your bam file is paired-end and correctly formed.\n'.format(total))
    elif unpaired > 0:
        sys.stderr.write('WARNING: it appears that {0} of the {1} reads overlapping the target STR regions were unpaired and so no useful data could be obtained from them.\n'.format(unpaired, total))

    # Write results
    if outfile:
        outstream = open(outfile, 'w')
    else:
        outstream = sys.stdout

    out_header = '\t'.join(['STR_chr', 'STR_start', 'STR_stop', 'motif', 'reflen', 'count'])
    outstream.write(out_header + '\n')
    summed.to_csv(outstream, sep='\t', header=False, index=False)
    outstream.close()

if __name__ == '__main__':
    main()
