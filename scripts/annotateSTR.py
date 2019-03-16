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
    """Convert pandas.DataFrame to BedTool object
    """
    df_asString = ['\t'.join([str(y) for y in x]) for x in df.values.tolist()]
    bed = bt.BedTool(df_asString)
    return(bed)

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
        if len(field_list) == 0 or field == '.':
            continue # ignore empty fields
        elif len(field_list) == 2:
            name = field_list[0].strip('"')
            value = field_list[1].strip('"')
            if value == '.' or value == '':
                value = 'NA'
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

def make_colnames(prefix, n):
    """ Make colnames in the form prefix+int. For example:
        prefix0, prefix1 ... prefix(n-1)
    Args:
        prefix (str): prefix to add to start of each name
        n (int): number of names to create
    Returns:
        list of strings
    """
    return([prefix + str(i) for i in range(n)])

def bt_annotate_df(target_df, annotation_file, annotation_colnames = None,
    tmp_bed = None):
    """Takes a pandas data frame (target_df), and returns it with new columns
    from an annotation supported by bedtools (e.g. bed/gff/gtf).
    Args:
        target_df (pandas.DataFrame): must be bed-compatible
        annotation_file (str): path to a gff/gtf file
        annotation_colnames (list): names for annotation_file columns (must
            match number of columns in annotation file)
        tmp_bed: path for bedtools temporary file. Default: randomly generated
            filename in form tmp-XXXXXXXX.bed in current directory.
    Returns:
        pandas.DataFrame
    """
    target_bed = dataframe_to_bed(target_df)
    target_colnames = list(target_df)
    # temp file for bedtools to write to and pandas to read
    if not tmp_bed:
        tmp_bed = 'tmp-' + randomletters(8) + '.bed'
    target_bed.intersect(b=annotation_file, loj=True, wb=True).saveas(tmp_bed)
    annotated_df = pd.read_csv(tmp_bed, sep='\t', header=None)
    os.remove(tmp_bed) #delete temporary file
    # if column names not provided for annotation file, make some up
    if not annotation_colnames:
        n_annotation_colnames = len(annotated_df.columns) - len(target_colnames)
        annotation_colnames = make_colnames('annotation_', n_annotation_colnames)
    colnames = target_colnames + annotation_colnames
    annotated_df.columns = colnames
    # Replace '.' (missing values from bedtools) with NA
    annotated_df.replace('.','NA', inplace=True)
    return(annotated_df)

def annotate_gff(target_df, annotation_file, annotation_cols = None,
    split_attribute = True):
    """Takes a pandas data frame (target_df) and returns it with a new column
    of annotation from a gff/gtf file (annotation_file).
    Args:
        target_df (pandas.DataFrame): must be bed-compatible
        annotation_file (str): path to a gff/gtf file
        annotation_cols (list): columns from annotation_file to include in
            output (default = all)
        split_attribute (bool): Replace attribute column multiple columns by
            splitting it on ';'
    Returns:
        pandas.DataFrame
    """
    gff_colnames = ['seqname', 'source', 'feature', 'ann_start', 'ann_end', 'score',
                    'strand', 'frame', 'attribute']
    target_colnames = list(target_df)

    annotated_df = bt_annotate_df(target_df, annotation_file,
        annotation_colnames = gff_colnames)

    #XXX check if there are any duplicate colnames and throw an error (need to
    # do with is other functions, so write a generic function to do it?)
    if annotation_cols:
        colnames_to_keep = target_colnames + annotation_cols
        annotated_df = annotated_df[colnames_to_keep]

    attribute_df = annotated_df['attribute'].apply(split_anntotation_col)
    if split_attribute:
        annotated_df = pd.concat([annotated_df, attribute_df], axis = 1)
        annotated_df = annotated_df.drop('attribute', axis=1)

    return annotated_df

def annotate_bed(target_df, bed_file, bed_colnames = None):
    """Takes a pandas data frame (target_df) and returns it with new column(s)
    of annotation from a bed file. Keeps the fourth and subsequent bed columns.
    Args:
        target_df (pandas.DataFrame): must be bed-compatible
        bed_file (str): path to a bed file
        bed_colnames (list): names of all bed columns (must match number of
            columns in bed file)
    Returns:
        pandas.DataFrame
    """
    annotated_df = bt_annotate_df(target_df, bed_file)

    # Fix up column names and remove those not required
    minimal_bed_colnames = ['bed_chrom', 'bed_start', 'bed_end']
    prefix = 'bed_annotation_' #XXX replace with filename?
    target_colnames = list(target_df)
    if bed_colnames:
        bed_colnames[:len(minimal_bed_colnames)] = minimal_bed_colnames
        colnames =  target_colnames + bed_colnames
    else:
        bed_colnames = minimal_bed_colnames
        n_bed_cols = len(annotated_df.columns) - len(target_colnames)
        bed_colnames = minimal_bed_colnames + make_colnames(prefix,
            n_bed_cols - len(minimal_bed_colnames))
        colnames = target_colnames + bed_colnames

    annotated_df.columns = colnames
    annotated_df = annotated_df.drop(minimal_bed_colnames, axis=1)

    return annotated_df

def dedup_annotations(str_df):
    """Remove duplicate rows from an annotated STR data frame, using the
    following rules:
    - CDS > 5' UTR > 3'UTR > exon > intron > other
    - take the first gene alphabetically
    - take the first transcript alphabetically
    - take the first pathogenic locus alphabetically
    Args:
        str_df (pandas.DataFrame): annotated STRetch results
    Returns:
        pandas.DataFrame
    """
    id_columns = ['chrom', 'start', 'end', 'sample', 'repeatunit']
    priority = ['CDS', 'start_codon', 'stop_codon', 'exon', 'intron']
    column_order = str_df.columns.values

    str_df['feature'] = pd.Categorical(str_df['feature'], priority)
    str_df.sort_values(by=['feature','gene_id','transcript_id', 'pathogenic'], inplace=True)
    str_df['feature'] = str_df['feature'].astype(str)
    str_df_dedup = str_df.groupby(id_columns).first().reset_index()
    return str_df_dedup[column_order]

def sortSTRs(str_df):
    """sort by outlier score then estimated size (bpInsertion), both descending
    Args:
        str_df (pandas.DataFrame): annotated STRetch results
    Returns:
        pandas.DataFrame
    """
    return(str_df.sort_values(['outlier', 'bpInsertion'], ascending=[False, False]))

def annotateSTRs(strfile, annfile, path_bed):
    """Take a STRetch results file and annotate it with a gene annotation file
    and pathogenic loci from a bed file.
    Args:
        strfile (str): path to a STRetch results file
        annfile (str): path to a gff/gtf file of gene annotations
        path_bed (str): path to a bed file containing pathogenic loci
    Returns:
        pandas.DataFrame
    """
    with open(strfile) as str_fhandle:
        str_df = pd.read_csv(str_fhandle, sep='\t')
        str_annotated = annotate_gff(str_df,
            annotation_file=annfile,
            annotation_cols=['feature','attribute'])

        if path_bed:
            str_annotated = annotate_bed(str_annotated, path_bed,
                bed_colnames=['bed_chrom', 'bed_start', 'bed_end', 'pathogenic'])
    #Dedup and sort
    str_annotated = dedup_annotations(str_annotated)
    str_annotated = sortSTRs(str_annotated)
    return str_annotated

def main(raw_args):
    # Parse command line arguments
    args = parse_args(raw_args)
    strfile = args.STRs_tsv
    pathfile = args.path_bed
    # genefile = args.genes
    annfile = args.annotation
    outfile = args.output

    annotated_df = annotateSTRs(strfile, annfile, pathfile)

    if outfile:
        annotated_df.to_csv(outfile, sep='\t', index = False)
    else:
        sys.stdout.write(annotated_df.to_csv(sep='\t', index = False))

if __name__ == '__main__':
    main(sys.argv[1:])
