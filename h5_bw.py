#!/usr/bin/env python
from optparse import OptionParser

import h5py
import numpy as np
import pyBigWig

'''
h5_bw.py

Convert a coverage HDF5 to BigWig.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <out_h5_file> <in_bw_file>'
    parser = OptionParser(usage)
    parser.add_option('-v', dest='verbose', default=False, action='store_true')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input HDF5 and output BigWig.')
    else:
        hdf5_file = args[0]
        bw_file = args[1]

    # open files
    h5_in = h5py.File(hdf5_file)
    bw_out = pyBigWig.open(bw_file, 'w')

    # construct header
    header = []
    chroms = sorted(h5_in.keys())
    for chrom in chroms:
        # chromosome and length
        header.append((chrom,len(h5_in[chrom])))

    # write header
    bw_out.addHeader(header)

    for chrom, length in header:
        if options.verbose:
            print(chrom)

        # read values
        x = np.array(h5_in[chrom])

        # write gzipped into HDF5
        bw_out.addEntries(chrom, 0, values=x, span=1, step=1)

    # close files
    h5_in.close()
    bw_out.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
