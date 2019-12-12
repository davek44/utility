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
    print_h5_tree(h5_in)

    # h5_keys = sorted(list(h5_in.keys()))
    # for hkey in h5_keys:
    #    print(h5_in[hkey])

    h5_in.close()

def print_h5_tree(h5_obj, depth=0):
    h5_keys = sorted(list(h5_obj.keys()))
    for hkey in h5_keys:
        print(''.join(['\t']*depth), h5_obj[hkey])
        if type(h5_obj[hkey]) == h5py._hl.group.Group:
            print_h5_tree(h5_obj[hkey], depth+1)
    
################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
