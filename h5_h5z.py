#!/usr/bin/env python
from optparse import OptionParser

import h5py
import numpy as np

from struct import pack, unpack

'''
h5_bw.py

Convert a coverage HDF5 to lossy compressed HDF5.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_h5_file> <out_h5_file>'
    parser = OptionParser(usage)
    parser.add_option('-v', dest='verbose', default=False, action='store_true')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input HDF5 and output BigWig.')
    else:
        in_h5_file = args[0]
        out_h5_file = args[1]

    # open files
    h5_in = h5py.File(in_h5_file)
    h5_out = h5py.File(out_h5_file, 'w')

    # construct header
    header = []
    chroms = sorted(h5_in.keys())
    for chrom in chroms:
        # chromosome and length
        header.append((chrom,len(h5_in[chrom])))

    for chrom, length in header:
        if options.verbose:
            print(chrom)

        # read values
        x = np.array(h5_in[chrom], dtype='float16')

        # write gzipped into HDF5
        h5_out.create_dataset(chrom, data=x, chunks=True, compression=32013, compression_opts=None, shuffle=False)

    # close files
    h5_in.close()
    h5_out.close()

def zfp_rate_opts(rate):
    """Create compression options for ZFP in fixed-rate mode

    The float rate parameter is the number of compressed bits per value.
    """
    ZFP_MODE_RATE = 1
    rate = pack('<d', rate)            # Pack as IEEE 754 double
    high = unpack('<I', rate[0:4])[0]  # Unpack high bits as unsigned int
    low = unpack('<I', rate[4:8])[0]   # Unpack low bits as unsigned int
    return (ZFP_MODE_RATE, 0, high, low, 0, 0)

def zfp_accuracy_opts(accuracy):
    """Create compression options for ZFP in fixed-accuracy mode

    The float accuracy parameter is the absolute error tolarance (e.g. 0.001).
    """
    ZFP_MODE_ACCURACY = 3
    accuracy = pack('<d', accuracy)        # Pack as IEEE 754 double
    high = unpack('<I', accuracy[0:4])[0]  # Unpack high bits as unsigned int
    low = unpack('<I', accuracy[4:8])[0]   # Unpack low bits as unsigned int
    print(high, low)
    return (ZFP_MODE_ACCURACY, 0, high, low, 0, 0)

def zfp_precision_opts(precision):
    """Create a compression options for ZFP in fixed-precision mode

    The float precision parameter is the number of uncompressed bits per value.
    """
    ZFP_MODE_PRECISION = 2
    return (ZFP_MODE_PRECISION, 0, precision, 0, 0, 0)

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
