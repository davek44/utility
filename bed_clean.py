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
    usage = 'usage: %prog [options] <bed_file> <sizes_file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != :
        parser.error('Must provide BED file and chrom sizes file')
    else:
        bed_file = args[0]
        csizes_file = args[1]

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
            shift = end - chrom_sizes[chrom]
            start -= shift
            end -= shift

            a[1] = str(start)
            a[2] = str(end)
            print '\t'.join(a)

    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
