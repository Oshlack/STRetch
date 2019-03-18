import sys
sys.path.append("..")
from annotateSTR import *
import pytest

# Define some input/output files that are shared between tests
str_file = 'test_data/test.tsv'
str_file_annotated = 'test_data.annotated.tsv'
annotation_file = 'test_data/gencode.v19.annotation.introns.gff3.gz'
annotation_mini = 'test_data/gencode.v19.annotation.introns.mini.gff3'
disease_bed = 'test_data/hg19.STR_disease_loci.bed'

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

def test_bt_annotate_df():
    target_df = pd.read_csv(str_file, sep='\t')
    annotated_df = bt_annotate_df(target_df, annotation_mini)

def test_annotate_gff():
    target_df = pd.read_csv(str_file, sep='\t')
    annotate_gff(target_df,
    annotation_mini,
    annotation_cols=['attribute'])

def test_annotate_bed():
    target_df = pd.read_csv(str_file, sep='\t')
    annotated_df = annotate_bed(target_df,
    bed_file=disease_bed,
    bed_colnames=['chrom', 'start', 'end', 'pathogenic']
    )

def test_annotateSTRs():
    bed_file=disease_bed
    str_annotated = annotateSTRs(str_file, annotation_mini, bed_file)
    str_annotated.to_csv(str_file_annotated, sep='\t', index = False)

def test_dedup_annotations():
    str_df = pd.read_csv(str_file_annotated, sep='\t')
    str_df_dedup = dedup_annotations(str_df)
    str_df_dedup.to_csv('tmp-dedup.tsv', sep='\t', index = False)
    assert len(str_df_dedup.index) == 2

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
