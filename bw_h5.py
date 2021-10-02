#!/usr/bin/env python
from optparse import OptionParser

import h5py
import numpy as np
import pyBigWig
import scipy.interpolate

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
    parser.add_option('-c', '--chr_strip', dest='chr_strip',
        default=False, action='store_true')
    parser.add_option('-i', dest='interp_nan',
        default=False, action='store_true',
        help='Interpolate NaNs [Default: %default]') 
    parser.add_option('-v', dest='verbose',
        default=False, action='store_true')
    parser.add_option('-z', dest='clip_zero',
        default=False, action='store_true',
        help='Clip negative values at zero [Default: %default]')
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
        x = bw_in.values(chrom, 0, chrom_lengths[chrom], numpy=True)

        # interpolate NaN
        if options.interp_nan:
            x = interp_nan(x)

        # clip negative values
        if options.clip_zero:
            x = np.clip(x, 0, np.inf)

        # strip "chr"
        if options.chr_strip:
            chrom = chrom.replace('chr','')

        # write gzipped into HDF5
        x = x.astype('float16')
        h5_out.create_dataset(chrom, data=x, dtype='float16', compression='gzip', shuffle=True)

    # close files
    h5_out.close()
    bw_in.close()


def interp_nan(x, kind='linear'):
    '''Linearly interpolate to fill NaN.'''

    # pad zeroes
    xp = np.zeros(len(x)+2)
    xp[1:-1] = x

    # find NaN
    x_nan = np.isnan(xp)

    if np.sum(x_nan) == 0:
        # unnecessary
        return x

    else:
        # interpolate
        inds = np.arange(len(xp))
        interpolator = scipy.interpolate.interp1d(
            inds[~x_nan],
            xp[~x_nan],
            kind=kind,
            bounds_error=False)

        loc = np.where(x_nan)
        xp[loc] = interpolator(loc)

        # slice off pad
        return xp[1:-1]

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
