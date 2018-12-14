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
    # Nonsense input - should probably generate errors XXX
    ((2,-10), (1,10), False),
])
def test_spans_region(test_region, target_region, expected):
    assert spans_region(test_region, target_region) == expected

def test_spans_region_errors(): #XXX incomplete
    test_region = (1,10)
    target_region = (4,6)
    assert spans_region(test_region, target_region)

def test_indel_size_nonspanning():
    """Raise ValueError if read doesn't span region"""
    bamfile = '11_L001_R1.STRdecoy.bam'
    test_read_name = '1-3871'
    region = (70713514, 70713561)
    chrom = 'chr13'
    bam = pysam.Samfile(bamfile, 'rb')
    with pytest.raises(ValueError):
        for read in bam.fetch():
            if read.query_name == test_read_name:
                indel_size(read, region, chrom)
                break

# def test_indel_size_0():
#     bamfile = '49_L001_R1.STRdecoy.bam'
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
#     bamfile = '49_L001_R1.STRdecoy.bam'
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
