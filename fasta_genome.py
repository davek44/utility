#!/usr/bin/env python
from optparse import OptionParser

import pysam

'''
Name

Description...
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fasta> <genome>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input FASTA and output genome files.')
    else:
        fasta_file = args[0]
        genome_file = args[1]

    fasta_open = pysam.Fastafile(fasta_file)
    genome_open = open(genome_file, 'w')
    
    for ref in fasta_open.references:
        ref_len = fasta_open.get_reference_length(ref)
        print('%s\t%d' % (ref, ref_len), file=genome_open)

    fasta_open.close()
    genome_open.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
