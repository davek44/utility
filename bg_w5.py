#!/usr/bin/env python
from optparse import OptionParser
import gzip

import h5py
import numpy as np

'''
bg_w5.py

Convert a BedGraph w/o overlapping entries to wig5.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_bg_file> <genome_file> <out_h5_file>'
    parser = OptionParser(usage)
    parser.add_option('-l', dest='norm_len',
        default=False, action='store_true',
        help='Normalize values by site length [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 3:
        parser.error('Must provide input BigWig, genome file, output HDF5.')
    else:
        bg_file = args[0]
        genome_file = args[1]
        hdf5_file = args[2]

    # initialize chromosome arrays
    chrm_values = {}
    for line in open(genome_file):
        a = line.split()
        chrm = a[0]
        chrm_len = int(a[1])
        chrm_values[chrm] = np.zeros(chrm_len, dtype='float16')

    # write bedgraph entries
    if bg_file[-3:] == '.gz':
        bg_open = gzip.open(bg_file, 'rt')
    else:
        bg_open = open(bg_file)

    for line in bg_open:
        if not line.startswith('#'):
            a = line.split()
            if len(a) >= 4:
                chrm = a[0]
                start = int(a[1])
                end = int(a[2])
                v = float(a[3])
                if options.norm_len:
                    v /= (end-start)
                chrm_values[chrm][start:end] = v

    bg_open.close()

    # write gzipped into HDF5
    h5_out = h5py.File(hdf5_file, 'w')
    for chrm in chrm_values:
        h5_out.create_dataset(chrm, data=np.nan_to_num(chrm_values[chrm]),
            dtype='float16', compression='gzip', shuffle=True)
    h5_out.close()



################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
