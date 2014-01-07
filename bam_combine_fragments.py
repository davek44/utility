#!/usr/bin/env python
from optparse import OptionParser
import pdb, sys
import pysam

################################################################################
# bam_combine_fragments.py
#
# Combine properly paired fragments in a BAM file to form a single entry.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()
    
    if len(args) == 1:
        bam_file = args[0]
    else:
        parser.error('Must provide BAM file')

    bam_in = pysam.Samfile(bam_file, 'rb')
    chromosomes = bam_in.references
    bam_out = pysam.Samfile('fragments.bam', 'wb', template=bam_in)
    bam_in.close()

    total_pairs = 0
    missing_pairs = 0
    overlapping = 0
    nonoverlapping = 0

    #for chrom_tid in range(len(chromosomes)):
    for chrom_tid in [1]:
        bam_in = pysam.Samfile(bam_file, 'rb')

        # hash properly paired reads by name to get read1 and read2 together
        proper_pairs = {}
        for aligned_read in bam_in:
            if aligned_read.tid == chrom_tid and aligned_read.mapq > 0 and aligned_read.is_proper_pair and aligned_read.rnext == aligned_read.tid:
                if not aligned_read.qname in proper_pairs:
                    proper_pairs[aligned_read.qname] = [{},{}]

                if aligned_read.is_read1:
                    read_num = 0
                else:
                    read_num = 1

                proper_pairs[aligned_read.qname][read_num][aligned_read.pos] = aligned_read

        # for each pair, walk the cigar strings simultaneously until they meet
        for qname in proper_pairs:
            for pos1, read1 in proper_pairs[qname][0].items():
                total_pairs += 1
                try:
                    read2 = proper_pairs[qname][1][read1.pnext]
                except:
                    missing_pairs += 1

                if read1.pos < read2.pos:
                    left_read = read1
                    right_read = read2
                else:
                    left_read = read2
                    right_read = read1

                if left_read.alen < right_read.pos:
                    nonoverlapping += 1
                else:
                    overlapping += 1

        bam_in.close()

    print >> sys.stderr, '%d of %d pairs missing.' % (missing_pairs,total_pairs)
    print >> sys.stderr, '%d overlapping, %d nonoverlapping' % (overlapping,nonoverlapping)
            
    
################################################################################
def tmp():
    for (operation,length) in aligned_read.cigar:
        # match
        if operation in [0,7,8]:
            if read_walked + length >= read_half:
                midpoint = genome_pos + (read_half - read_walked)
                break
            else:
                genome_pos += length
                read_walked += length

        # insertion
        elif operation in [1,3]:
            genome_pos += length

        # deletion
        elif operation == 2:
            read_walked += length


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
