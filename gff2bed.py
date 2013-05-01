#!/usr/bin/env python
from optparse import OptionParser
import gff

################################################################################
# gff2bed.py
#
# Convert a gff file to a bed file. Each entry is converted independently,
# so no blocks.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gff file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide gff file')
    else:
        gff_file = args[0]

    for line in open(gff_file):
        a = line.split('\t')

        cols = [a[0], str(int(a[3])-1), a[4], a[2], '0', a[6], '0', '0', '255,0,0', '1', str(int(a[4])-int(a[3])+1), '0']
        print '\t'.join(cols)
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
