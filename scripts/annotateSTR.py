#!/usr/bin/env python
import argparse
import random
import string
import pybedtools as bt
import pandas as pd
import sys
import os

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = argparse.ArgumentParser(description='Annotate an STRetch results file (STRs.tsv) with gene/transcript information')
    parser.add_argument(
        'STRs_tsv', type=str,
        help='Input STRs.tsv file produced by STRetch.')
    parser.add_argument(
        '--path', type=str,
        help='BED file containing the locations of known pathogenic STR loci.')
    parser.add_argument(
        '--genes', type=str,
        help='BED file containing the locations of genes.')
    parser.add_argument(
        '--annotations', type=str, nargs="+",
        help='One or more BED files containing annotations such as exons, introns etc.')
    parser.add_argument(
        '--output', type=str, required=False,
        help='Output file name. Defaults to stdout.')

    return parser.parse_args()

def randomletters(length):
   return ''.join(random.choice(string.ascii_lowercase) for i in range(length))

def annotate_bed(target_df, annotation_bed, column_heading=None):
    """Takes a pandas data frame (target_df) of STRetch output, and returns it 
        with a new column of annotation from the bed file (annotation_bed).
        '.' if no annotation available at that locus. Only the first column 
        in the bed file (after the positions) is used as annotation."""

    if not column_heading:
        column_heading = annotation_bed
    #convert data frame to list of strings so that bedtools can read it in
    target_df_asString = ['\t'.join([str(y) for y in x]) for x in target_df.values.tolist()]
    target_bed = bt.BedTool(target_df_asString)
    target_colnames = list(target_df)

    with open(annotation_bed) as annotation_fhandle:
        annotation_bed = bt.BedTool(annotation_fhandle)
        tmp_bed = 'tmp-' + randomletters(8) + '.bed' #temp file for bedtools to write to and pandas to read
        target_bed.intersect(b=annotation_bed, loj=True).saveas(tmp_bed)
        colnames = target_colnames + ['match_chr', 'match_start', 'match_end', column_heading]
        annotated_df = pd.read_csv(tmp_bed, sep='\t', header=None)
        annotated_df = annotated_df.iloc[:,:len(colnames)] # remove following columns
        annotated_df.columns = colnames
        annotated_df.drop(['match_chr', 'match_start', 'match_end'], axis=1, inplace=True)
        os.remove(tmp_bed) #delete temporary file

    return annotated_df


def main():
    # Parse command line arguments
    args = parse_args()
    strfile = args.STRs_tsv
    pathfile = args.path
    genefile = args.genes
    annfiles = args.annotations
    outfile = args.output

    if outfile:
        outstream = open(outfile, 'w')
    else:
        outstream = sys.stdout

    # Create an annotation master bed file
    # Use bed arithmetic to get a file with a single annotation (exon, intron 
    # etc) for each position in the genome
    
    ann_names = [x.split('.bed.sorted')[0].split('_')[-1] for x in annfiles]
    ann_beds = [bt.BedTool(bedfile) for bedfile in annfiles]

#    bedfile1 = annfiles[0]
#    ann_name1 = bedfile1.split('.bed.sorted')[0].split('_')[-1]
#    bed1 = bt.BedTool(bedfile1)
#
#    bedfile2 = annfiles[1]
#    ann_name2 = bedfile2.split('.bed.sorted')[0].split('_')[-1]
#    bed2 = bt.BedTool(bedfile2)

    # Subtract bed files from each other 
    cat_bed = bt.BedTool.cat(*ann_beds[:-1], postmerge=True, force_truncate=True)
    ann_beds[0].merge().subtract(cat_bed)
    sys.exit()
    bed_sub = ann_beds[0].merge().subtract(*bt.BedTool.cat(ann_beds[:-1], 
                                    postmerge=True, force_truncate=True))
    bed_sub.saveas('combined_annotation.bed')
    #bed2_sub = ann_beds[1].merge().subtract(ann_beds[0].merge())
    #bed2_sub.saveas('combined_annotation.bed')

#    with open(strfile) as str_fhandle:
#        str_df = pd.read_csv(str_fhandle, sep='\t')
#   
#        str_df_annotated = annotate_bed(str_df, pathfile, "pathlocus")
#        str_df_annotated = annotate_bed(str_df_annotated, genefile, "gene")
#
#        str_df_annotated.to_csv(outfile, sep='\t')

if __name__ == '__main__':
    main()
