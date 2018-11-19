#!/usr/bin/env python
from optparse import OptionParser
import pdb

import numpy as np
import pysam

'''
make_fasta_genome.py

Make a "genome" file, with chromosome names and lengths.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fasta_file> <genome_file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input FASTA file and output genome file')
    else:
        fasta_file = args[0]
        genome_file = args[1]

    genome_out = open(genome_file, 'w')

    fasta_open = pysam.Fastafile(fasta_file)
    ref_indexes = np.argsort(fasta_open.lengths)[::-1]

    for i in ref_indexes:
        print('%s\t%d' % (fasta_open.references[i], fasta_open.lengths[i]), file=genome_out)
    fasta_open.close()

    genome_out.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
