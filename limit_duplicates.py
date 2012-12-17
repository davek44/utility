#!/usr/bin/env python
from optparse import OptionParser
import pdb, sys
import pysam

################################################################################
# limit_duplicates.py
#
# Accept a BAM file as input and max out the number of reads that can occur
# at one single position.
#
# The script is imperfect in terms of handling a PCR duplicate that maps to
# many positions in the genome, but if the alignments are sorted, saving the
# removed reads should minimize the issue.
#
# Lately, I've preferred samtools rmdup which only considers chrom, strand
# and start point.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam file>'
    parser = OptionParser(usage)
    parser.add_option('-d', dest='dup_t', type='int', default=10, help='Number of duplicates at which removal begins [Default: %default]')
    parser.add_option('-o', dest='output_bam', default='out.bam', help='Ouput BAM file')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide bam file')
    else:
        bam_file = args[0]

    chr_aligns = {}
    removed_reads = set()

    bam_in = pysam.Samfile(bam_file, 'rb')
    bam_out = pysam.Samfile(options.output_bam, 'wb', header=bam_in.header)

    for read in bam_in:
        rkey = (read.tid,read.pos,read.alen)
        if not read.qname in removed_reads:
            if chr_aligns.get(rkey,0) < options.dup_t:
                bam_out.write(read)
                chr_aligns[rkey] = chr_aligns.get(rkey,0) + 1
            else:
                removed_reads.add(read.qname)

    bam_in.close()
    bam_out.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
