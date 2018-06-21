#!/usr/bin/env python
from optparse import OptionParser

import h5py
import numpy as np
import zarr

'''
h5_zarr.py

Convert a coverage HDF5 to BigWig.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_h5_file> <out_zarr_file>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='chunk_size', default=None, type='int')
    parser.add_option('-v', dest='verbose', default=False, action='store_true')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input HDF5 and output BigWig.')
    else:
        hdf5_file = args[0]
        zarr_file = args[1]

    # open files
    h5_in = h5py.File(hdf5_file)
    zarr_out = zarr.open_group(zarr_file, 'w')

    # foreach chromosome
    for chrom in h5_in.keys():
        if options.verbose:
            print(chrom)

        # read values
        x = np.array(h5_in[chrom], dtype='float16')

        # write gzipped into HDF5
        z = zarr_out.create_dataset(chrom, data=x, shape=x.shape, dtype='float16', chunks=options.chunk_size)
        if options.verbose:
            print(z)

    # close files
    h5_in.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
