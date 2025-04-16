#!/usr/bin/env python
from optparse import OptionParser
import sys

import h5py
import numpy as np
import pyBigWig
import scipy.interpolate

'''
bw_w5.py

Convert a BigWig to wigHDF5.
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
    parser.add_option('-m', dest='min_norm',
        default=False, action='store_true',
        help='Normalize the minimum nonzero value to 1 [Default: %default]')
    parser.add_option('-s', dest='scale',
        default=1.0, type='float',
        help='Scale all values (e.g. to undo normalization) [Default: %default]')
    parser.add_option('-v', dest='verbose',
        default=False, action='store_true')
    parser.add_option('-z', dest='clip_zero',
        default=False, action='store_true',
        help='Clip negative values at zero [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input BigWig and output HDF5.')
    else:
        bw_files = args[0]
        hdf5_file = args[1]

    # open files
    bw_files = bw_files.split(',')
    bw_ins = [pyBigWig.open(bw_file) for bw_file in bw_files]
    h5_out = h5py.File(hdf5_file, 'w')

    # process chromosomes in length order
    chrom_lengths = bw_ins[0].chroms()
    chroms = sorted(chrom_lengths.keys())
    length_chroms = [(chrom_lengths[chrm],chrm) for chrm in chroms]
    length_chroms = sorted(length_chroms)[::-1]
    min_factor = None

    # for each chromosome
    for clength, chrom in length_chroms:
        if options.verbose:
            print(chrom)

        # read values
        x = bw_ins[0].values(chrom, 0, clength, numpy=True)
        for bw_in in bw_ins[1:]:
            x += bw_in.values(chrom, 0, clength, numpy=True)

        # scale
        if options.scale != 1:
            x = x*options.scale

        # normalize min to 1
        #  (a simple strategy to undo normalization)
        if options.min_norm:
            if min_factor is None:
                min_factor = x[x>0].min()
                print('Min normalization factor: %f' % min_factor, file=sys.stderr)
            x /= min_factor

        # interpolate NaN
        if options.interp_nan:
            x = interp_nan(x)
        else:
            x = np.nan_to_num(x)

        # clip negative values
        if options.clip_zero:
            x = np.clip(x, 0, np.inf)

        # clip float16 min/max
        x = np.clip(x, np.finfo(np.float16).min, np.finfo(np.float16).max)
        x = x.astype('float16')

        # strip "chr"
        if options.chr_strip:
            chrom = chrom.replace('chr','')

        # write gzipped into HDF5
        h5_out.create_dataset(chrom, data=x, dtype='float16', compression='gzip', shuffle=True)

    # close files
    h5_out.close()
    for bw_in in bw_ins:
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
