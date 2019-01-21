import sys
sys.path.append("..")
from decoy_STR import *
import pytest

@pytest.mark.parametrize("sequence, expected", [
    ('A', ['A']),
    ('AT', ['AT', 'TA']),
    ('CAG', ['CAG', 'AGC', 'GCA']),
])
def test_circular_permuted(sequence, expected):
    assert circular_permuted(sequence) == expected

@pytest.mark.parametrize("sequence, expected", [
    ('A', ['A', 'T']),
    ('AT', ['AT', 'AT']),
    ('CAG', ['CAG', 'CTG']),
])
def test_self_and_rev_complement(sequence, expected):
    assert self_and_rev_complement(sequence) == expected

@pytest.mark.parametrize("sequence, expected", [
    ('T', 'A'),
    ('CAG', 'AGC'),
    ('AAT', 'AAT'),
])
def test_normalise_str(sequence, expected):
    assert normalise_str(sequence) == expected

def test_randomletters_len():
    outstream = open('test_decoySTR.output.txt', 'w')
    write_decoys(outstream, seqlength = 20, repeatunits = None)

def test_write_decoys_seqlength():
    f = 'test_decoySTR.output.txt'
    with open(f, 'w') as outstream:
        seqlength = 10
        repeatunits = None
        write_decoys(outstream, seqlength, repeatunits)
    with open(f, 'r') as instream:
        all_lines = instream.readlines()
        assert '>STR-A\n' in all_lines
        assert len(all_lines) == 501 * 2
        assert len(all_lines[1].strip()) == seqlength

def test_write_decoys_repeatunits():
    f = 'test_decoySTR.output.txt'
    with open(f, 'w') as outstream:
        seqlength = None
        repeatunits = 5
        write_decoys(outstream, seqlength, repeatunits)
    with open(f, 'r') as instream:
        all_lines = instream.readlines()
        assert '>STR-A\n' in all_lines
        assert len(all_lines) == 501 * 2
        assert len(all_lines[5].strip()) == repeatunits * 2

def test_parse_args_defaults():
    """Check correct defaults are set when no options/arguments given"""
    args = parse_args([])
    assert args.output == None
    assert args.length == 500
    assert args.repeatunits == None
