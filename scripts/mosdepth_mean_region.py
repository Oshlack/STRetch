#!/usr/bin/env python
"""Estimate the mean coverage from the mosdepth region.bed.gz coverage output
"""

import argparse
import sys
import pandas as pd
import numpy as np
import gzip

__author__ = "Harriet Dashnow"
__credits__ = ["Harriet Dashnow"]
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "h.dashnow@gmail.com"

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description='Estimate the mean coverage from the mosdepth region.bed.gz coverage output.')
    parser.add_argument(
        'REGION_BED', type=str, 
        help='region.bed.gz mosdepth output file (gzipped bed file mean coverages for each region)')
    parser.add_argument(
        '--out', type=str, required=False,
        help='Output file name. Defaults to stdout.')
    return parser.parse_args()

def parse_bed(filename):
    """Parse mosdepth bed with coverage to pandas df"""
    try:
        data = pd.read_table(filename, delim_whitespace = True)
        if data.shape[1] == 4:
            data.columns = ['chr', 'start', 'stop', 'coverage']
        elif data.shape[1] == 5:
            data.columns = ['chr', 'start', 'stop', 'region', 'coverage']
        else:
            sys.exit('ERROR: file {0} had an unexpected number of columns (not 3 or 4).\n'.format(filename))
    except pd.io.common.EmptyDataError:
        sys.exit('ERROR: file {0} was empty.\n'.format(filename))
    return(data)

def main():
    args = parse_args()
    region_bed = args.REGION_BED
    outfile = args.out

    if outfile:
        outstream = open(outfile, 'w')
    else:
        outstream = sys.stdout

    # Parse input file
    region_data = parse_bed(region_bed)

    # Calculate mean
    #region_median = np.median(region_data['coverage'])
    region_mean = np.mean(region_data['coverage'])

    outstream.write(str(region_mean)+'\n')

if __name__ == '__main__':
    main()
