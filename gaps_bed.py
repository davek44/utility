#!/usr/bin/env python
from optparse import OptionParser
import pdb

################################################################################
# gaps_bed.py
#
# Print a bed file of the gaps in a fasta file.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fasta file>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='gap_size', default=50, type='int', help='Minimum gap size to print a bed entry [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide fasta file')
    else:
        fasta_file = args[0]

    for line in open(fasta_file):
        if line[0] == '>':
            chrom = line[1:].rstrip()
            seq_i = 0
            gap_start = None
        else:
            for nt in line.rstrip():
                if nt == 'N':
                    if gap_start == None:
                        gap_start = seq_i                        
                else:
                    if gap_start != None and seq_i-gap_start >= options.gap_size:
                        print '\t'.join([chrom,str(gap_start),str(seq_i)])
                    gap_start = None

                seq_i += 1


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
