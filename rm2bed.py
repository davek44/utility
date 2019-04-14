#!/usr/bin/env python
from optparse import OptionParser
import gzip

'''
rm2bed.py

Convert RepeatMasker .out format to BED.
'''

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
            rm_in = gzip.open(args[0], 'rt')
        else:
            rm_in = open(args[0])

    for i in range(4):
        line = rm_in.readline()
    while line:
        a = line.split()

        chrm = a[4]
        start = str(int(a[5])-1)
        end = a[6]

        if a[8] == '+':
            strand = '+'
        else:
            strand = '-'

        repeat = a[9]
        family = a[10]

        cols = (chrm, start, end, '%s;%s' % (family,repeat), '.', strand)
        print('\t'.join(cols))

        line = rm_in.readline()

    rm_in.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
