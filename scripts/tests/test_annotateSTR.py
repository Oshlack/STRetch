import sys
sys.path.append("..")
from annotateSTR import *
import pytest

# Currently working with hg19_gencodeV19_comp.gtf
# other formats not working
@pytest.mark.parametrize("annotation, expected", [
    ('gene_id ""ENST00000414504.2""; transcript_id ""ENST00000414504.2""; ',
        {'gene_id': 'ENST00000414504.2', 'transcript_id': 'ENST00000414504.2'} ),
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
    target_df = pd.read_csv('test_data/test.tsv', sep='\t')
    annotation_file = 'test_data/hg19_gencodeV19_comp.gtf.gz'
    annotated_df = bt_annotate_df(target_df, annotation_file)

def test_annotate_gff():
    target_df = pd.read_csv('test_data/test.tsv', sep='\t')
    annotate_gff(target_df,
    annotation_file='test_data/hg19_gencodeV19_comp.gtf.gz',
    #tmp_bed='tmp.bed',
    annotation_cols=['attribute'])

def test_annotate_bed():
    target_df = pd.read_csv('test_data/test.tsv', sep='\t')
    annotated_df = annotate_bed(target_df,
    bed_file='test_data/hg19.STR_disease_loci.bed',
    bed_colnames=['chrom', 'start', 'end', 'pathogenic']
    )

def test_annotateSTRs():
    str_file = 'test_data/test.tsv'
    annotation_file = 'test_data/hg19_gencodeV19_comp.gtf.gz'
    bed_file='test_data/hg19.STR_disease_loci.bed'
    annotateSTRs(str_file, annotation_file, bed_file)
