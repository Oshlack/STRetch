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

def test_annotate_gff():
    data_dir = '/Users/harriet.dashnow/Documents/git/STR-variation/data/'
    target_df = pd.read_csv(data_dir + 'test.tsv', sep='\t')
    annotate_gff(target_df,
    annotation_file=data_dir + 'hg19_gencodeV19_comp.gtf.gz',
    tmp_bed='tmp.bed',
    annotation_cols=['attribute'])

def test_annotateSTRs():
    data_dir = '/Users/harriet.dashnow/Documents/git/STR-variation/data/'
    str_file = data_dir + 'test.tsv'
    annotation_file = data_dir + 'hg19_gencodeV19_comp.gtf.gz'
    annotateSTRs(str_file, annotation_file)
