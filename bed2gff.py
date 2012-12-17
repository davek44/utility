#!/usr/bin/env python
from optparse import OptionParser
import gff

################################################################################
# bed2gff.py
#
# Convert a bed file to a gff file. No blocks.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bed file>'
    parser = OptionParser(usage)
    parser.add_option('--source', dest='source', default='bed', help='Gff format "source" [Default: %default]')
    parser.add_option('--feature', dest='feature', default='feature', help='Gff format "feature" [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide bed file')
    else:
        bed_file = args[0]

    group_num = 0
    for line in open(bed_file):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        if len(a) >= 5:
            score = a[4]
        else:
            score = '.'
        if len(a) >= 6:
            strand = a[5]
        else:
            strand = '+'
        group_num += 1

        cols = [a[0], options.source, options.feature, str(int(a[1])+1), a[2], score, strand, '.', str(group_num)]
        print '\t'.join(cols)        
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
