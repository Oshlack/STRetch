#!/usr/bin/env python
import argparse
import random
import string
import pandas as pd
import sys
import os
import gzip

import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import pybedtools as bt


def parse_args(raw_args):
    "Parse the input arguments, use '-h' for help"
    parser = argparse.ArgumentParser(description='Annotate an STRetch results file (STRs.tsv) with gene/transcript information')
    parser.add_argument(
        'STRs_tsv', type=str,
        help='Input STRs.tsv file produced by STRetch.')
    parser.add_argument(
        '--path_bed', type=str,
        help='BED file containing the locations of known pathogenic STR loci.')
    # parser.add_argument(
    #     '--genes', type=str,
    #     help='BED file containing the locations of genes.')
    parser.add_argument(
        '--annotation', type=str,
        help='GTF/GFF genome annotation file (can be gzipped)')
    parser.add_argument(
            '--ids', type=str,
            help='Gene IDs and their corresponding gene names/symbols') #XXX describe format
    parser.add_argument(
        '--output', type=str,
        help='Output file name. Defaults to stdout.')

    return parser.parse_args(raw_args)

def randomletters(length):
   return ''.join(random.choice(string.ascii_lowercase) for i in range(length))

def dataframe_to_bed(df):
    """Convert pandas dataframe to BedTool object
    """
    df_asString = ['\t'.join([str(y) for y in x]) for x in df.values.tolist()]
    bed = bt.BedTool(df_asString)
    return(bed)

def annotate_bed(target_df, annotation_bed, column_heading=None):
    """Takes a pandas data frame (target_df) of STRetch output, and return it
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

# Currently working with hg19_gencodeV19_comp.gtf
# other formats not working
def parse_gff_annotation(annotation):
    """
    Currently expecting only two columns in the annotation
    Args:
        annotation (str): the string of a single line in the GTF annotation
            column
    Returns:
        pandas.DataFrame
    """
    annotation_dict = {}
    fields = annotation.split(';')
    for field in fields:
        # Strip whitespace
        field_list = [x.strip() for x in field.split() if not x == '']
        if len(field_list) == 0:
            continue # ignore empty fields
        elif len(field_list) == 2:
            name = field_list[0].strip('"')
            value = field_list[1].strip('"')
            annotation_dict[name] = value
        else:
            raise Exception(('unexpected number of attributes in field '
                                '(expecting 2): {}').format(field))

    return annotation_dict

def split_anntotation_col(annotation):
    ann_dict = parse_gff_annotation(annotation)
    #ann_df = pd.DataFrame(ann_dict)
    #return(ann_df)
    #return([ann_dict[i][0] for i in ann_dict])
    return(pd.Series(ann_dict))

def bt_annotate_df(target_df, annotation_file):
    """Takes a pandas data frame (target_df), and returns it with new columns
    from an annotation supported by bedtools (e.g. bed/gff/gtf).
    Args:
        target_df (pandas.DataFrame): must be bed-compatible
        annotation_file (str): path to a gff/gtf file
        annotation_cols (list): columns from annotation_file to include in
            output (default = all)
        split_attribute (bool): Replace attribute column multiple columns by
            splitting it on ';'
        null: null value to use, will be '.' if no annotation available at that locus
        tmp_bed: path for bedtools temporary file. Default: randomly generated
            filename in form tmp-XXXXXXXX.bed in current directory.
    Returns:
        pandas.DataFrame
    """
    pass


def annotate_gff(target_df, annotation_file, annotation_cols = None,
    split_attribute = True, null = '.', tmp_bed = None):
    """Takes a pandas data frame (target_df) and returns it with a new column
    of annotation from a gff/gtf file (annotation_file).
    Args:
        target_df (pandas dataframe): must be bed-compatible
        annotation_file (str): path to a gff/gtf file
        annotation_cols (list): columns from annotation_file to include in
            output (default = all)
        split_attribute (bool): Replace attribute column multiple columns by
            splitting it on ';'
        null: null value to use, will be '.' if no annotation available at that locus
        tmp_bed: path for bedtools temporary file. Default: randomly generated
            filename in form tmp-XXXXXXXX.bed in current directory.
    Returns:
        pandas.DataFrame
    """
    gff_colnames = ['seqname', 'source', 'feature', 'start', 'end', 'score',
                    'strand', 'frame', 'attribute']

    target_bed = dataframe_to_bed(target_df)
    target_colnames = list(target_df)

    # temp file for bedtools to write to and pandas to read
    if not tmp_bed:
        tmp_bed = 'tmp-' + randomletters(8) + '.bed'
    target_bed.intersect(b=annotation_file, loj=True, wb=True).saveas(tmp_bed)
    colnames = target_colnames + gff_colnames
    annotated_df = pd.read_csv(tmp_bed, sep='\t', header=None)
    annotated_df.columns = colnames
    os.remove(tmp_bed) #delete temporary file

    if annotation_cols:
        colnames_to_keep = target_colnames + annotation_cols
        annotated_df = annotated_df[colnames_to_keep]

    attribute_df = annotated_df['attribute'].apply(split_anntotation_col)
    if split_attribute:
        annotated_df = pd.concat([annotated_df, attribute_df], axis = 1)
        annotated_df = annotated_df.drop('attribute', axis=1)

    return annotated_df

#annotated_df.to_csv('tmp.tsv', sep='\t', index=False)

def annotateSTRs(strfile, annfile):

    with open(strfile) as str_fhandle:
        str_df = pd.read_csv(str_fhandle, sep='\t')
        str_annotated = annotate_gff(str_df,
            annotation_file=annfile,
            annotation_cols=['attribute'])
    return str_annotated



    # Create an annotation master bed file
    # Use bed arithmetic to get a file with a single annotation (exon, intron
    # etc) for each position in the genome

    # ann_names = [x.split('.bed.sorted')[0].split('_')[-1] for x in annfiles]
    # ann_beds = [bt.BedTool(bedfile) for bedfile in annfiles]

#    bedfile1 = annfiles[0]
#    ann_name1 = bedfile1.split('.bed.sorted')[0].split('_')[-1]
#    bed1 = bt.BedTool(bedfile1)
#
#    bedfile2 = annfiles[1]
#    ann_name2 = bedfile2.split('.bed.sorted')[0].split('_')[-1]
#    bed2 = bt.BedTool(bedfile2)

    # Subtract bed files from each other
    # cat_bed = bt.BedTool.cat(*ann_beds[:-1], postmerge=True, force_truncate=True)
    # ann_beds[0].merge().subtract(cat_bed)
    # sys.exit()
    # bed_sub = ann_beds[0].merge().subtract(*bt.BedTool.cat(ann_beds[:-1],
    #                                 postmerge=True, force_truncate=True))
    # bed_sub.saveas('combined_annotation.bed')
    #bed2_sub = ann_beds[1].merge().subtract(ann_beds[0].merge())
    #bed2_sub.saveas('combined_annotation.bed')

#    with open(strfile) as str_fhandle:
#        str_df = pd.read_csv(str_fhandle, sep='\t')
#
#        str_df_annotated = annotate_bed(str_df, pathfile, "pathlocus")
#        str_df_annotated = annotate_bed(str_df_annotated, genefile, "gene")
#
#        str_df_annotated.to_csv(outfile, sep='\t')

def main(raw_args):
    # Parse command line arguments
    args = parse_args(raw_args)
    strfile = args.STRs_tsv
    # pathfile = args.path
    # genefile = args.genes
    annfile = args.annotation
    outfile = args.output

    if outfile:
        outstream = open(outfile, 'w')
    else:
        outstream = sys.stdout

    outstream.write(annotateSTRs(strfile, annfile))
    outstream.close()



if __name__ == '__main__':
    main(sys.argv[1:])
