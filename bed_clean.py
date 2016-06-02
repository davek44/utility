#!/usr/bin/env python
from optparse import OptionParser

################################################################################
# bed_clean.py
#
# Detect and correct BED regions extending beyond chromosome ends
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <csizes_file> <bed_file>'
    parser = OptionParser(usage)
    parser.add_option('-d', dest='delete', default=False, action='store_true', help='Delete entries beyond boundaries [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide BED file and chrom sizes file')
    else:
        csizes_file = args[0]
        bed_file = args[1]

    # read in chromosome sizes
    chrom_sizes = {}
    for line in open(csizes_file):
        a = line.split()
        chrom_sizes[a[0]] = int(a[1])

    # clean BED file
    for line in open(bed_file):
        a = line.split()
        chrom = a[0]
        start = int(a[1])
        end = int(a[2])

        if end <= chrom_sizes[chrom]:
            print line,

        else:
            if not options.delete:
                end = chrom_sizes[chrom]
                a[2] = str(end)
                print '\t'.join(a)
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
