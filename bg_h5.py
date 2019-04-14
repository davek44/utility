#!/usr/bin/env python
from optparse import OptionParser

import h5py
import numpy as np

'''
bg_h5.py

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
    for line in open(genome_file):
        a = line.split()
        chrm = a[0]
        chrm_len = int(a[1])
        chrm_values[chrm] = np.zeros(chrm_len, dtype='float16')

    # write bedgraph entries
    for line in open(bg_file):
        a = line.split()
        chrm = a[0]
        start = int(a[1])
        end = int(a[2])
        v = float(a[3])
        chrm_values[chrm][start:end] += v

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
