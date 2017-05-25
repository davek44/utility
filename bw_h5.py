#!/usr/bin/env python
from optparse import OptionParser

import h5py
import numpy as np
import pyBigWig

'''
bw_h5.py

Convert a BigWig to HDF5.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_bw_file> <out_h5_file>'
    parser = OptionParser(usage)
    parser.add_option('-v', dest='verbose', default=False, action='store_true')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input BigWig and output HDF5.')
    else:
        bw_file = args[0]
        hdf5_file = args[1]

    # open files
    bw_in = pyBigWig.open(bw_file)
    h5_out = h5py.File(hdf5_file, 'w')

    # for each chromosome
    chrom_lengths = bw_in.chroms()
    for chrom in chrom_lengths:
        if options.verbose:
            print(chrom)

        # read values
        x = bw_in.values(chrom, 0, chrom_lengths[chrom], numpy=True).astype('float16')

        # write gzipped into HDF5
        h5_out.create_dataset(chrom, data=x, dtype='float16', compression='gzip', shuffle=True)

    # close files
    h5_out.close()
    bw_in.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
