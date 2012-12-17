#!/usr/bin/env python
from optparse import OptionParser
import glob, sys, pdb

################################################################################
# nuc2gff.py
#
# Convert alignments from a nucmer coords file to features in a .gff file.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <coords file>'
    parser = OptionParser(usage)
    parser.add_option('-f', dest='feature_name', default='domain', help='Feature name [Default: %default]')
    parser.add_option('-p', dest='pct_t', type='float', default=0.9, help='Percentage of the 2nd sequence that must be covered by the alignment [Default: %default]')
    parser.add_option('-i', dest='idy_t', type='float', default=0.8, help='% identity that must be exceeded by the alignment [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide nucmer output coords file')
    else:
        coords_file = args[0]

    # get header
    cf = open(coords_file)
    for i in range(5):
        cf.readline()

    line = cf.readline()
    while line:
        a = line.split()
        
        header1 = a[-2]
        header2 = a[-1]

        start1 = int(a[0])
        end1 = int(a[1])

        idy = float(a[9])/100.0
        len2 = int(a[12])
        if int(a[3]) < int(a[4]):
            strand = '+'
            start2 = int(a[3])
            end2 = int(a[4])
        else:
            strand = '-'
            start2 = int(a[4])
            end2 = int(a[3])

        if idy > options.idy_t and end2-start2+1 > len2*options.pct_t:
            gff_a = [header1, 'nucmer', options.feature_name, str(start1), str(end1), '.', strand, '.', header2]
            print '\t'.join(gff_a)
        
        line = cf.readline()
 

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
