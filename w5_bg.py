#!/usr/bin/env python
from optparse import OptionParser

import h5py
import numpy as np

'''
w5_bg.py

Convert a Wiggle HDF5 to BedGraph.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_w5> <out_bg>'
    parser = OptionParser(usage)
    parser.add_option('-v', dest='verbose', default=False, action='store_true')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input Wig5 and output BedGraph.')
    else:
        in_w5_file = args[0]
        out_bg_file = args[1]

    # open files
    in_w5_open = h5py.File(in_w5_file)
    out_bg_open = open(out_bg_file, 'w')

    header = 'track type=bedGraph'
    print(header, file=out_bg_open)

    for chrm in sorted(in_w5_open.keys()):
        if options.verbose:
            print(chrm, flush=True)

        # read values
        x = np.array(in_w5_open[chrm])

        # write to bedgraph
        i = 0
        while i < len(x):
            start = i
            end = i+1
            while end < len(x) and x[start] == x[end]:
                end += 1

            cols = [chrm, str(start), str(end), '%.4f'%x[start]]
            print('\t'.join(cols), file=out_bg_open)

            i = end

    in_w5_open.close()
    out_bg_open.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
