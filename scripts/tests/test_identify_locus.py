import sys
sys.path.append("..")
from identify_locus import *
import pytest
import pysam

def test_randomletters_len():
    length = 5
    letters = randomletters(length)
    assert len(letters) == length

def test_randomletters_diff():
    length = 10
    letters1 = randomletters(length)
    letters2 = randomletters(length)
    assert letters1 != letters2

# Genotypes
# ==> 11 <==
# chr13 70713514 70713561 CTG_16/201
# ==> 49 <==
# chr13 70713514 70713561 CTG_16/1

#fetch('chr1', 100, 120)

@pytest.mark.parametrize("test_region, target_region, expected", [
    # Easy ones
    ((1,10), (4,6), True),
    ((1,2), (4,6), False),
    # Borderline cases
    ((1,10), (1,10), True),
    ((1,9), (1,10), False),
    ((2,10), (1,10), False),
    # Single base
    ((1,10), (5,5), True),
    ((1,10), (10,10), True),
    ((1,10), (1,1), True),
    # Nonsense input - should probably generate errors XXX
    ((2,-10), (1,10), False),
])
def test_spans_region(test_region, target_region, expected):
    assert spans_region(test_region, target_region) == expected

@pytest.mark.parametrize("test_region, target_region", [
    # Nonsense input - should generate errors
    ((2,-10), (1)),
    ((2,-10), (1,3,6)),
])
def test_spans_region_errors(test_region, target_region):
    with pytest.raises(TypeError):
        spans_region(test_region, target_region)

@pytest.mark.parametrize("position, region, expected", [
    # Easy ones
    (5, (1,10), True),
    (20, (1,10), False),
    # Borderline cases
    (1, (1,10), True),
    (10, (1,10), True),
    (0, (1,10), False),
    (11, (1,10), False),
    # Nonsense input - should probably generate errors XXX
    #((1,2), (1,10), False),
])
def test_in_region(position, region, expected):
    assert in_region(position, region) == expected

@pytest.mark.parametrize("position, region", [
    # Nonsense input - should generate errors
    (2, (1,10,20)),
])
def test_in_region_errors(position, region):
    with pytest.raises(TypeError):
        in_region(position, region)

def test_indel_size_nonspanning():
    """Raise ValueError if read doesn't span region"""
    bamfile = 'test_data/11_L001_R1.STRdecoy.bam'
    test_read_name = '1-3871'
    region = (70713514, 70713561)
    chrom = 'chr13'
    bam = pysam.Samfile(bamfile, 'rb')
    with pytest.raises(ValueError):
        for read in bam.fetch():
            if read.query_name == test_read_name:
                indel_size(read, region, chrom)
                break

@pytest.mark.parametrize("test_read_name, expected", [
    ('1-293', 0), # same as ref
    ('1-293:10I', 10), # easy insertion
    ('1-293:5D', -5), # easy deletion
    ('1-293:0compound', 0), # insertion and deletion same size cancel out
    ('1-293:2Dcompound', -2), # insertion and deletion of different sizes makes deletion
    ('1-293:20Icompound', 20), # two insertions added together
    ('1-293:outside', 0), # indel outside STR
])
def test_indel_size(test_read_name, expected):
    """Test indel corretly identified"""
    bamfile = 'test_data/49_tests.STRdecoy.sam'
    region = (70713514, 70713561)
    chrom = 'chr13'
    bam = pysam.Samfile(bamfile, 'rb')
    for read in bam.fetch():
        if read.query_name == test_read_name:
            try:
                assert indel_size(read, region, chrom) == expected
            except ValueError:
                continue

# def test_indel_size_0():
#     bamfile = 'test_data/49_L001_R1.STRdecoy.bam'
#     test_read_name = '1-293'
#     region = (70713514, 70713561)
#     fetch_region = (70713514,)
#     chrom = 'chr13'
#     bam = pysam.Samfile(bamfile, 'rb')
#     i = 0
#     for read in bam.fetch(chrom, *fetch_region):
#         #print(read.cigar)
#         if read.query_name == test_read_name:
#             this_indel_size = indel_size(read, region, chrom)
#             print(read.cigar)
#             #break
#             print(read)
# 
#     #assert False
#     #assert this_indel_size == 0
# 
# 
# 
# def test_indel_size():
#     bamfile = 'test_data/49_L001_R1.STRdecoy.bam'
#     test_read_name = '1-293'
#     region = (70713514, 70713561)
#     fetch_region = (70713514,)
#     chrom = 'chr13'
#     bam = pysam.Samfile(bamfile, 'rb')
#     i = 0
#     for read in bam.fetch(chrom, *fetch_region):
#         #print(read.cigar)
#         try: 
#             print(indel_size(read, region, chrom))
#             print(read.cigar)
#             #break
#             if read.cigar[0][1] == 150:
#                 print(read)
#             #i += 1
#             
#             #if i > 1000:
#             #    break
#         except ValueError:
#             continue
# #    assert False
