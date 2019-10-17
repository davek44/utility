#!/usr/bin/env python
from optparse import OptionParser
import gzip

import h5py
import numpy as np

'''
bg_w5.py

Convert a BedGraph to HDF5.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_bg_file> <genome_file> <out_h5_file>'
    parser = OptionParser(usage)
    parser.add_option('-v', dest='verbose', default=False, action='store_true')
    (options,args) = parser.parse_args()

    if len(args) != 3:
        parser.error('Must provide input BigWig, genome file, output HDF5.')
    else:
        bg_file = args[0]
        genome_file = args[1]
        hdf5_file = args[2]

    # initialize chromosome arrays
    chrm_values = {}
    chrm_counts = {}
    for line in open(genome_file):
        a = line.split()
        chrm = a[0]
        chrm_len = int(a[1])
        chrm_values[chrm] = np.zeros(chrm_len, dtype='float16')
        chrm_counts[chrm] = np.zeros(chrm_len, dtype='uint8')

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
                chrm_values[chrm][start:end] += v
                chrm_counts[chrm][start:end] += 1

    bg_open.close()

    # take mean
    for chrm in chrm_values:
        chrm_values[chrm] = np.divide(chrm_values[chrm], chrm_counts[chrm])
        chrm_values[chrm] = np.nan_to_num(chrm_values[chrm])

    # write gzipped into HDF5
    h5_out = h5py.File(hdf5_file, 'w')
    for chrm in chrm_values:
        h5_out.create_dataset(chrm, data=chrm_values[chrm], dtype='float16', compression='gzip', shuffle=True)
    h5_out.close()



################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
