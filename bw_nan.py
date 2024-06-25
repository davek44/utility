#!/usr/bin/env python
from optparse import OptionParser

import numpy as np
import pyBigWig

'''
bw_nan.py

Compute NaN % in a BigWig.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_bw_file> <out_h5_file>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()
    if len(args) != 1:
        parser.error('Must provide input BigWig.')
    else:
        bw_file = args[0]

    # open files
    bw_in = pyBigWig.open(bw_file)

    # process chromosomes in length order
    chrom_lengths = bw_in.chroms()
    chroms = sorted(chrom_lengths.keys())
    length_chroms = [(chrom_lengths[chrm],chrm) for chrm in chroms]
    length_chroms = sorted(length_chroms)[::-1]
    mode_factor = None

    total_nt = 0
    nan_nt = 0

    # for each chromosome
    for clength, chrom in length_chroms:
        # read values
        x = bw_in.values(chrom, 0, chrom_lengths[chrom], numpy=True)

        # find NaN
        x_nan = np.isnan(x)

        total_nt += len(x)
        nan_nt += x_nan.sum()

    # close files
    bw_in.close()

    nan_pct = nan_nt / total_nt
    print('%.6f' % nan_pct)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
