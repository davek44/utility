#!/usr/bin/env python
from optparse import OptionParser

import h5py
import numpy as np
import zarr

'''
zarr_h5.py

Convert a coverage Zarr to HDF5.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_zarr_file> <out_h5_file>'
    parser = OptionParser(usage)
    parser.add_option('-v', dest='verbose', default=False, action='store_true')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input Zarr and output HDF5.')
    else:
        zarr_file = args[0]
        hdf5_file = args[1]

    # open files
    zarr_in = zarr.open_group(zarr_file, 'r')
    h5_out = h5py.File(hdf5_file, 'w')

    # foreach chromosome
    for chrom in sorted(zarr_in.keys()):
        if options.verbose:
            print(chrom)

        # read values
        x = np.array(zarr_in[chrom])

        # write gzipped into HDF5
        h5_out.create_dataset(chrom, data=x, dtype='float16', chunks=True, compression='lzf', shuffle=True)

    # close files
    h5_out.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
