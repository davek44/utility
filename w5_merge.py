#!/usr/bin/env python
from optparse import OptionParser

import h5py
import numpy as np

'''
w5_merge.py

Merge wig5 files using a specified summary statistic.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <out_w5> <in1_w5> <in2_w5> ...'
    parser = OptionParser(usage)
    parser.add_option('-s', dest='sum_stat',
        default='sum', help='Summary statistic [Default: %default]')
    parser.add_option('-v', dest='verbose',
        default=False, action='store_true')
    (options,args) = parser.parse_args()

    if len(args) < 3:
        parser.error('Must provide output and two or more input wig5.')
    else:
        out_w5_file = args[0]
        in_w5_files = args[1:]

    # open input wig5
    in_w5_opens = [h5py.File(iwf) for iwf in in_w5_files]
    in_num = len(in_w5_opens)

    # take keys union
    in_keys = set()
    for in_w5_open in in_w5_opens:
        in_keys |= in_w5_open.keys()

    # open output file
    out_w5_open = h5py.File(out_w5_file, 'w')

    for out_key in in_keys:
        if options.verbose:
            print(out_key)

        # read data
        in_key_len = len(in_w5_opens[0][out_key])
        in_key_data = np.zeros((in_num,in_key_len), dtype='float32')
        for i in range(in_num):
            in_key_data[i] = np.array(in_w5_opens[i][out_key])

        # summarize
        if options.sum_stat == 'sum':
            out_key_data = in_key_data.sum(axis=0)
        elif options.sum_stat == 'mean':
            out_key_data = in_key_data.mean(axis=0)
        else:
            print('Cannot identify summary statistic %s' % options.sum_stat)

        # write
        out_w5_open.create_dataset(out_key, data=out_key_data.astype('float16'), dtype='float16')

    out_w5_open.close()



################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
