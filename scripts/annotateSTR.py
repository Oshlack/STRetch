#!/usr/bin/env python
import argparse
import random
import string
import pandas as pd
import sys
import os
import gzip
import toolz

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
    parser.add_argument(
        '--omim', type=str,
        help='OMIM mim2gene.txt file (can be downloaded from https://www.omim.org/downloads/)')
    parser.add_argument(
        '--annotation', type=str,
        help='GFF3 genome annotation file (can be gzipped)')
    parser.add_argument(
        '--tss', type=str,
        help='GFF3 genome annotation file of transcription start sites (TSS). Will be caculated from annotation if not provided (can be gzipped)')
    parser.add_argument(
        '--output', type=str,
        help='Output file name. Defaults to stdout.')
    parser.add_argument(
        '--chunksize', type=int, default = 10000,
        help='Annotate this many variants at a time to reduce memory requirements')

    return parser.parse_args(raw_args)

def randomletters(length):
   return ''.join(random.choice(string.ascii_lowercase) for i in range(length))

@toolz.curry
def remove_prefix(s, prefix):
    """if string s starts with prefix, return s without the prefix,
    otherwise return original s
    """
    if s.startswith(prefix):
        return(s[len(prefix):])
    else:
        return s

@toolz.curry
def add_prefix(s, prefix):
    """if string s doesn't start with prefix, return s with the prefix,
    otherwise return original s
    """
    if s.startswith(prefix):
        s
    else:
        return prefix + s

def dataframe_to_bed(df):
    """Convert pandas.DataFrame to BedTool object
    """
    df_asString = ['\t'.join([str(y) for y in x]) for x in df.values.tolist()]
    bed = bt.BedTool(df_asString)
    return(bed)

def parse_gff_annotation(annotation):
    """
    Parse the annotation column of a GFF3 format genome annotation file into a
    pandas data frame
    Args:
        annotation (str): annotation column from a single line in a GFF3 file
    Returns:
        pandas.DataFrame
    """
    annotation_dict = {}
    fields = annotation.split(';')
    for field in fields:
        # Strip whitespace
        field_list = [x.strip() for x in field.split('=') if not x == '']
        if len(field_list) == 0 or field == '.' or field == 'NA':
            continue # ignore empty fields
        elif len(field_list) == 1:
            sys.stderr.write('Unexpectedly short field: ' + field + '\n')
            continue # ignore single fields
        else:
            name = field_list[0].strip('"')
            value = field_list[1].strip('"')
            if value == '.' or value == '':
                value = 'NA'
            annotation_dict[name] = value
    return annotation_dict

def split_anntotation_col(annotation):
    ann_dict = parse_gff_annotation(annotation)
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

def calculate_TSS(gff_line):
    """Calculate the position of the transcription start site (TSS) from a
    single transcript annotation line from a gff3 file
    Args:
        gff (str): transcript line from a gff3 file
    Returns:
        str: a new gff3 format line describing the TSS
    """
    splitline = gff_line.split()
    if not splitline[2] == 'transcript':
        pass #XXX raise error?
    start = splitline[3]
    end = splitline[4]
    strand = splitline[6]
    #attributes = splitline[8]
    if strand == '+':
        tss = start
    elif strand == '-':
        tss = end
    else:
        tss = None # or 'NA'?
    #attributes_dict = parse_gff_annotation(attributes)
    #transcript_id = attributes_dict['transcript_id']
    #gene_id = attributes_dict['gene_id']
    #gene_name = attributes_dict['gene_name']
    return_list = splitline
    return_list[1] = 'annotateSTR'
    return_list[2] = 'TSS'
    return_list[3] = return_list[4] = tss
    return_string = '\t'.join(return_list)
    return return_string

def open_check_gz(filename):
    """If filename ends in .gz open with gzip, else open normally.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def gff_TSS(gff_in, gff_out, gff_tmp = None):
    """Calculate the positions of transcription start sites (TSSs) from a gff3
    file and write them to a bed file.
    Args:
        gff_in (str): path to a gff3 input file containing transcripts
        gff_out (str): path to a gff3 output file (will be overwritten)
        gff_tmp: path for bedtools gff temporary file. Default: randomly generated
            filename in form tmp-XXXXXXXX.gff in the ./tmp directory.
    """
    if not gff_tmp:
        tmpdir = 'tmp'
        gff_tmp = tmpdir + '/tmp-' + randomletters(8) + '.gff'
        try:
            os.mkdir(tmpdir)
        except FileExistsError:
            pass

    with open_check_gz(gff_in) as fin, open(gff_tmp, 'w') as fout:
        gff_header = """##gff-version 3
#description: TSSs calculated as the first nucleotide in each transcript
#format: gff3
"""
        fout.write(gff_header)
        for line in fin:
            if not line.startswith('#'):
                splitline = line.split()
                if splitline[2] == 'transcript':
                    fout.write(calculate_TSS(line)+'\n')

        # bedtools sort -i test_tss.gff -header
        gff = bt.BedTool(gff_tmp)
        gff.sort(header=True).saveas(gff_out)
        os.remove(gff_tmp) #delete temporary file

def bt_annotate_df(target_df, annotation_file, command = 'intersect',
    annotation_colnames = None, tmp_bed = None):
    """Takes a pandas data frame (target_df), and returns it with new columns
    from an annotation from a gff3 file.
    Args:
        target_df (pandas.DataFrame): must be bed-compatible
        annotation_file (str): path to a gff3 file
        command (str): BedTools command to use - 'intersect' or 'closest'
        annotation_colnames (list): names for annotation_file columns (must
            match number of columns in annotation file)
        tmp_bed: path for bedtools temporary file. Default: randomly generated
            filename in form tmp-XXXXXXXX.bed in the ./tmp directory.
    Returns:
        pandas.DataFrame
    """
    target_bed = dataframe_to_bed(target_df)
    target_colnames = list(target_df)
    # temp file for bedtools to write to and pandas to read
    if not tmp_bed:
        tmpdir = 'tmp'
        tmp_bed = tmpdir + '/tmp-' + randomletters(8) + '.bed'
        try:
            os.mkdir(tmpdir)
        except FileExistsError:
            pass

    if command == 'intersect':
        target_bed.intersect(b=annotation_file, loj=True, wb=True).saveas(tmp_bed)
    elif command == 'closest':
        target_bed.sort().closest(b=annotation_file, D='b', t='first').saveas(tmp_bed)
    else:
        raise ValueError('Unknown command: ' + command +
            ', valid options are "intersect" or "closest"')

    annotated_df = pd.read_csv(tmp_bed, sep='\t', header=None)
    os.remove(tmp_bed) #delete temporary file

    #XXX check for the situation that there were no matches, causes missing columns

    if annotation_colnames:
        duplicate_colnames = [colname for colname in annotation_colnames
            if colname in target_colnames]
        if len(duplicate_colnames) > 0:
            raise ValueError('Column names must be unique, these were duplicated: '
                + str(duplicate_colnames))
    # if column names not provided for annotation file, make some up
    else:
        n_annotation_colnames = len(annotated_df.columns) - len(target_colnames)
        annotation_colnames = make_colnames('annotation_', n_annotation_colnames)

    colnames = target_colnames + annotation_colnames
    annotated_df.columns = colnames
    # Replace '.' (missing values from bedtools) with NA
    annotated_df.replace('.','NA', inplace=True)
    return(annotated_df)

def annotate_gff(target_df, annotation_file, command = 'intersect',
    keep_cols = None, split_attributes = True, col_prefix = 'new_'):
    """Takes a pandas data frame (target_df) and returns it with a new column(s)
    of annotation from a gff3 file (annotation_file).
    Args:
        target_df (pandas.DataFrame): must be bed-compatible
        annotation_file (str): path to a gff3 file
        command (str): BedTools command to use - 'intersect' or 'closest'
        keep_cols (list): columns from annotation_file to include in
            output. Can include named values from the gff attributes columns
            if 'attributes' is included, then all named values from the attributes
            column are kept (default is None which results in keeping all cols)
        split_attributes (bool): Replace attribute column with multiple columns
            by splitting it on ';'
        col_prefix (str): prefix to use to make column names unique if required
    Returns:
        pandas.DataFrame
    """
    gff_colnames = ['seqname', 'source', 'feature', 'ann_start', 'ann_end', 'score',
                    'strand', 'frame', 'attributes']
    attributes_colname = 'attributes'
    if command == 'closest':
        gff_colnames.append('distance')
    target_colnames = list(target_df)

    # check for and rename duplicate colnames
    duplicate_colnames = [colname for colname in gff_colnames
        if colname in target_colnames]

    # Add a prefix to the gff column names to make them unique
    gff_colnames = [col_prefix + colname for colname in gff_colnames]
    attributes_colname = col_prefix + attributes_colname
    keep_cols = [col_prefix + colname for colname in keep_cols]

    annotated_df = bt_annotate_df(target_df, annotation_file, command,
        annotation_colnames = gff_colnames)

    attribute_df = annotated_df[attributes_colname].apply(split_anntotation_col)
    # Rename attribute columns with prefix
    add_col_prefix = add_prefix(prefix=col_prefix)
    attribute_df = attribute_df.rename(add_col_prefix, axis='columns')

    if split_attributes:
        annotated_df = pd.concat([annotated_df, attribute_df], axis = 1)
        annotated_df = annotated_df.drop(attributes_colname, axis=1)

    if keep_cols:
        if attributes_colname in keep_cols:
            keep_cols += list(attribute_df)
            keep_cols.remove(attributes_colname)
        colnames_to_keep = target_colnames + keep_cols
        annotated_df = annotated_df[colnames_to_keep]

    remove_col_prefix = remove_prefix(prefix=col_prefix)
    annotated_df = annotated_df.rename(remove_col_prefix, axis='columns')

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
    priority = ['CDS', 'start_codon', 'stop_codon',
        'stop_codon_redefined_as_selenocysteine', 'exon', 'UTR', 'intron',
        'transcript', 'gene']
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

def annotate_tss(str_df, tss_gff):
    """
    """
    tss_annotated = annotate_gff(str_df, tss_gff, command = 'closest',
        keep_cols=['transcript_name', 'distance'], col_prefix = 'tss_')
    # Rename TSS columns
    tss_annotated = tss_annotated.rename(
        columns={'transcript_name': 'closest_TSS_transcript_name',
        'distance': 'closest_TSS_distance'})
    return tss_annotated

def parse_omim(omim_file):
    """Parse mim2gene.txt file downloaded from https://www.omim.org/downloads/
    Extracts out genes only
    Args:
        omim_file (str): path to OMIM file, e.g. mim2gene.txt
    Returns:
        pandas.DataFrame
    """
    omim_df = pd.read_csv(omim_file, sep='\t', comment='#',
        names = ['mim_number', 'mim_entry_type', 'entrez_gene_id',
        'hgnc_gene_symbol', 'ensembl_gene_id'])
    return omim_df.loc[omim_df['mim_entry_type'] == 'gene']

#XXX make path_bed optional?
def annotateSTRs(str_df, annfile, path_bed, tss_file, omim_file=None):
    """Take a STRetch results file and annotate it with a gene annotation file
    and pathogenic loci from a bed file.
    Args:
        str_df (pandas.DataFrame): annotated STRetch results
        annfile (str): path to a gff3 file of gene annotations
        path_bed (str): path to a bed file containing pathogenic loci
        tss_file (str): path to a gff3 file of TSS positions
        omim_file (str): path to OMIM file, e.g. mim2gene.txt
    Returns:
        pandas.DataFrame
    """
    str_annotated = annotate_gff(str_df,
        annotation_file=annfile,
        keep_cols=['feature','gene_name','gene_id','transcript_id'], col_prefix='gff_')
    str_annotated.to_csv('another-tmp.tsv', sep='\t', index = False)

    # Check if annotated genes are in OMIM
    if omim_file:
        omim_df = parse_omim(omim_file)
        # check if gene_id in omim_df['ensembl_gene_id']
        base_gene_id = str_annotated['gene_id'].apply(lambda x: str(x).split('.')[0])
        gene_id_in_omim = base_gene_id.isin(omim_df['ensembl_gene_id'])
        # check if gene_name in omim_df['hgnc_gene_symbol']
        gene_name_in_omim = str_annotated['gene_name'].isin(omim_df['hgnc_gene_symbol'])
        # either gene_id or gene_name in omim
        str_annotated['in_omim'] = gene_id_in_omim | gene_name_in_omim
        # check if TSS is in omim? if so move this to below tss annotation

    # annotate with TSSs
    str_annotated = annotate_tss(str_annotated, tss_gff=tss_file)

    # Annotate pathogenic loci
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
    omim_file = args.omim
    annfile = args.annotation
    outfile = args.output
    chunk_size = args.chunksize

    #XXX allow these to be provdied at commandline?
    tss_file = None

    # if no TSS file provided, calculate TSS positions from transcripts
    if not tss_file:
        tss_file = 'calculated_TSS.gff' #XXX better file name here (or require it be set?)
        gff_TSS(annfile, tss_file)

    with open(strfile) as str_fhandle:
        str_reader = pd.read_csv(str_fhandle, sep='\t', chunksize=chunk_size)
        header_written = False
        for str_df in str_reader:
            annotated_df = annotateSTRs(str_df, annfile, pathfile,
                tss_file, omim_file)

            # need to delete file if it exists (to avoid appending an old one)
            # and write header for first chunk only
            if not header_written:
                if outfile:
                    annotated_df.to_csv(outfile, sep='\t', index = False,
                        header = True, mode = 'w')
                else:
                    sys.stdout.write(annotated_df.to_csv(sep='\t', index = False,
                        header = True))
                header_written = True
            else:
                if outfile:
                    annotated_df.to_csv(outfile, sep='\t', index = False,
                        header = False, mode = 'a')
                else:
                    sys.stdout.write(annotated_df.to_csv(sep='\t', index = False,
                        header = False))

if __name__ == '__main__':
    main(sys.argv[1:])
