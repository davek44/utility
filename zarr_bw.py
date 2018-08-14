#!/usr/bin/env python
from optparse import OptionParser

import zarr
import numpy as np
import pyBigWig

'''
zarr_bw.py

Convert a coverage Zarr to BigWig.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_zarr_file> <out_bw_file>'
    parser = OptionParser(usage)
    parser.add_option('-v', dest='verbose', default=False, action='store_true')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input HDF5 and output BigWig.')
    else:
        zarr_file = args[0]
        bw_file = args[1]

    # open files
    zarr_in = zarr.open_group(zarr_file, 'r')
    bw_out = pyBigWig.open(bw_file, 'w')

    # construct header
    header = []
    chroms = sorted(zarr_in.keys())
    for chrom in chroms:
        # chromosome and length
        header.append((chrom,len(zarr_in[chrom])))

    # write header
    bw_out.addHeader(header)

    for chrom, length in header:
        if options.verbose:
            print(chrom)

        # read values
        x = np.array(zarr_in[chrom])

        # write gzipped into HDF5
        bw_out.addEntries(chrom, 0, values=x, span=1, step=1)

    # close files
    bw_out.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
