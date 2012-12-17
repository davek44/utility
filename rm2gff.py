#!/usr/bin/env python
from optparse import OptionParser
import gzip
import gff

################################################################################
# rm2gff.py
#
# Convert RepeatMasker .out format to gff
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <rm out>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide RepeatMasker .out file')
    else:
        if args[0][-2:] == 'gz':
            rm_in = gzip.open(args[0])
        else:
            rm_in = open(args[0])

    for i in range(4):
        line = rm_in.readline()
    while line:
        a = line.split()

        if a[8] == '+':
            strand = '+'
        else:
            strand = '-'

        cols = (a[4], 'RepeatMasker', 'repeat', a[5], a[6], '.', strand, '.', gff.kv_gtf({'repeat':a[9], 'family':a[10]}))
        print '\t'.join(cols)

        line = rm_in.readline()
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
