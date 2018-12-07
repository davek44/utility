#!/usr/bin/env python
from optparse import OptionParser
import h5py

'''
Name

Description...
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <h5_file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    h5_in = h5py.File(args[0], 'r')

    h5_keys = sorted(list(h5_in.keys()))
    for hkey in h5_keys:
        print(h5_in[hkey])

    h5_in.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
