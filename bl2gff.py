#!/usr/bin/env python
from optparse import OptionParser
import glob, sys, pdb

################################################################################
# bl2gff.py
#
# Convert alignments from my Blast output to features in a .gff file.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <blast file>'
    parser = OptionParser(usage)
    parser.add_option('-f', dest='feature_name', default='domain', help='Feature name [Default: %default]')
    parser.add_option('-p', dest='pct_t', type='float', default=0, help='Percentage of the 2nd sequence that must be covered by the alignment [Default: %default]')
    parser.add_option('-i', dest='idy_t', type='float', default=0, help='% identity that must be exceeded by the alignment [Default: %default]')    
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide blast output file')
        exit(1)
    else:
        blast_file = args[0]

    for line in open(blast_file):
        a = line.split()
        
        header1 = a[-2]
        header2 = a[-1]

        start1 = int(a[0])
        end1 = int(a[1])

        alen2 = int(a[7])
        len2 = int(a[10])
        idy = float(a[12])
        if int(a[3]) < int(a[4]):
            strand = '+'
        else:
            strand = '-'

        if idy > options.idy_t and alen2 > len2*options.pct_t:
            gff_a = [header1, 'blast', options.feature_name, str(start1), str(end1), '.', strand, '.', header2]
            print '\t'.join(gff_a)
        

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
