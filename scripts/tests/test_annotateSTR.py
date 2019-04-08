import sys
sys.path.append("..")
from annotateSTR import *
import pytest
import numpy as np

# Define some input/output files that are shared between tests
str_file = 'test_data/test.tsv'
str_file_annotated = 'test_data.annotated.tsv'
annotation_file = 'test_data/gencode.v19.annotation.introns.gff3.gz'
annotation_mini = 'test_data/gencode.v19.annotation.introns.mini.gff3'
tss_file = 'test_data/gencode.v19.annotation.TSS.gff'
disease_bed = 'test_data/hg19.STR_disease_loci.bed'
#omim_file = 'test_data/mim2gene.txt'
omim_file = 'test_data/biomart_GRCh37.p13_omim.tsv'

@pytest.mark.parametrize("s, prefix, expected", [
    ('new_colname', 'new_', 'colname'),
    ('colname', 'new_', 'colname'),
])
def test_remove_prefix(s, prefix, expected):
    assert remove_prefix(s, prefix) == expected

@pytest.mark.parametrize("annotation, expected", [
    ('gene_id=ENSG00000230223.2;transcript_id=ENST00000414504.2',
        {'gene_id': 'ENSG00000230223.2', 'transcript_id': 'ENST00000414504.2'} ),
    ('gene_id=.;transcript_id=.',
        {'gene_id': 'NA', 'transcript_id': 'NA'} ),
    # ('gene_id ""someval""; transcript_id val1 val2; ',
    #     {'gene_id': ['someval'], 'transcript_id': ['val1', 'val2']} ),
])
def test_parse_gff_annotation(annotation, expected):
    assert parse_gff_annotation(annotation) == expected

def test_parse_args_none():
    """Exit and display usage when no args/options are provided"""
    with pytest.raises(SystemExit):
        parser = parse_args([])

@pytest.mark.parametrize("prefix, n, expected", [
    ('col_', 2, ['col_0', 'col_1']),
])
def test_make_colnames(prefix, n, expected):
    assert make_colnames(prefix, n) == expected

@pytest.mark.parametrize("gff_line, expected", [
    ('chr1	HAVANA	transcript	11869	14409	.	+	.	ID=transcript1;Parent=gene1;gene_id=ENSG00000223972.4;transcript_id=ENST00000456328.2;gene_type=pseudogene;gene_status=KNOWN;gene_name=DDX11L1;transcript_type=processed_transcript;transcript_status=KNOWN;transcript_name=DDX11L1-002;level=2;havana_gene=OTTHUMG00000000961.2;havana_transcript=OTTHUMT00000362751.1;tag=basic',
        'chr1	annotateSTR	TSS	11869	11869	.	+	.	ID=transcript1;Parent=gene1;gene_id=ENSG00000223972.4;transcript_id=ENST00000456328.2;gene_type=pseudogene;gene_status=KNOWN;gene_name=DDX11L1;transcript_type=processed_transcript;transcript_status=KNOWN;transcript_name=DDX11L1-002;level=2;havana_gene=OTTHUMG00000000961.2;havana_transcript=OTTHUMT00000362751.1;tag=basic' ),
])
def test_calculate_TSS(gff_line, expected):
    assert calculate_TSS(gff_line) == expected

def test_gff_TSS():
    output_gff = 'test_tss.gff'
    gff_TSS(annotation_mini, output_gff)
    with open(output_gff) as f:
        for line in f:
            if not line.startswith('#'):
                assert line.startswith('chr1	annotateSTR	TSS	11869	11869')
                break

# def test_gff_TSS_all():
#     output_gff = 'transcripts_tss.gff'
#     gff_TSS(annotation_file, output_gff)

def test_bt_annotate_df():
    target_df = pd.read_csv(str_file, sep='\t')
    annotated_df = bt_annotate_df(target_df, annotation_mini)

def test_bt_annotate_df_closest():
    target_df = pd.read_csv(str_file, sep='\t')
    annotated_df = bt_annotate_df(target_df, 'test_tss.gff', command = 'closest')
    #annotated_df.to_csv('tmp-closest-tss.tsv', sep='\t', index = False)

def test_bt_annotate_df_dup():
    """Raise error if duplicate colnames"""
    target_df = pd.read_csv('test_data/test.annotated_w_duplicates.tsv', sep='\t')
    with pytest.raises(ValueError):
        annotated_df = bt_annotate_df(target_df, tss_file,
            annotation_colnames = ['seqname', 'source', 'feature', 'ann_start',
            'ann_end', 'score', 'strand', 'frame', 'attributes'])

def test_annotate_gff():
    target_df = pd.read_csv(str_file, sep='\t')
    annotated_df = annotate_gff(target_df,
    annotation_mini, keep_cols=['attributes'])
    #annotated_df.to_csv('tmp-annotate_gff.tsv', sep='\t', index = False)

def test_annotate_gff_closest():
    target_df = pd.read_csv(str_file, sep='\t')
    annotated_df = annotate_gff(target_df,
    tss_file, command = 'closest',
    keep_cols=['transcript_name', 'distance'])
    annotated_df.to_csv('tmp-closest-tss.tsv', sep='\t', index = False)

#XXX Test the situation that there's no matching results
# def test_annotate_gff_closest_none():
#     target_df = pd.read_csv(str_file, sep='\t')
#     annotated_df = annotate_gff(target_df,
#     'test_tss.gff', command = 'closest',
#     keep_cols=['transcript_name', 'distance'])
#     annotated_df.to_csv('tmp-closest-tss.tsv', sep='\t', index = False)

def test_annotate_tss():
    target_df = pd.read_csv(str_file, sep='\t')
    annotated_df = annotate_tss(target_df, tss_gff=tss_file)
    annotated_df.to_csv('tmp-annotate-tss.tsv', sep='\t', index = False)

def test_annotate_bed():
    target_df = pd.read_csv(str_file, sep='\t')
    annotated_df = annotate_bed(target_df,
    bed_file=disease_bed,
    bed_colnames=['chrom', 'start', 'end', 'pathogenic']
    )

def test_parse_omim():
    omim_df = parse_omim('test_data/mim2gene.txt')
    omim_df.to_csv('omim_genes.tsv', sep='\t', index = False)

def test_parse_biomart_omim():
    omim_df = parse_biomart_omim('test_data/biomart_GRCh37.p13_omim.tsv')
    #print(omim_df['ensembl_gene_id'][:2])
    assert np.array_equal(omim_df['ensembl_gene_id'][:3],
        ['ENSG00000261577', 'ENSG00000261577', 'ENSG00000261258'])

def test_in_omim():
    str_df = pd.read_csv(str_file, sep='\t')
    str_annotated = annotateSTRs(str_df, annotation_mini, disease_bed, tss_file = tss_file)
    str_annotated_omim = in_omim(str_annotated, omim_file)
    print(str_annotated)
    print(str_annotated_omim)
    str_annotated_omim['in_omim'] == [True, True, False]


def test_annotateSTRs():
    bed_file=disease_bed
    str_df = pd.read_csv(str_file, sep='\t')
    str_annotated = annotateSTRs(str_df, annotation_mini, bed_file, tss_file, omim_file)
    str_annotated.to_csv(str_file_annotated, sep='\t', index = False)

def test_dedup_annotations():
    str_df = pd.read_csv(str_file_annotated, sep='\t')
    str_df_dedup = dedup_annotations(str_df)
    str_df_dedup.to_csv('tmp-dedup.tsv', sep='\t', index = False)
    assert len(str_df_dedup.index) == 3

@pytest.mark.parametrize("str_dict, expected", [
    ({'outlier': [5, 10], 'bpInsertion': [100, 40]},
    {'outlier': [10, 5], 'bpInsertion': [40, 100]}),
    # sort by bpInsertion if outlier missing or identical
    ({'outlier': [None, None], 'bpInsertion': [10, 40]},
    {'outlier': [None, None], 'bpInsertion': [40, 10]}),
])
def test_sortSTRs(str_dict, expected):
    str_df = pd.DataFrame(data=str_dict)
    str_df_sorted = sortSTRs(str_df)
    str_df_sorted.to_csv('tmp-sorted.tsv', sep='\t', index = False)
    assert str_df_sorted.reset_index(drop=True).equals(pd.DataFrame(expected))
