#!/usr/bin/env python
from optparse import OptionParser
import pdb, os
import pysam

################################################################################
# bam_plus_minus.py
#
# Separate the alignments in a BAM file into two BAM files of the plus and
# minus strand reads.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam file>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide bam file')
    else:
        bam_file = args[0]

    chr_starts = {}

    bam_pre = os.path.splitext(bam_file)[0]

    bam_in = pysam.Samfile(bam_file, 'rb')
    bamp_out = pysam.Samfile('%s_p.bam'%bam_pre, 'wb', header=bam_in.header)
    bamm_out = pysam.Samfile('%s_m.bam'%bam_pre, 'wb', header=bam_in.header)

    for read in bam_in:
        if read.is_reverse:
            bamm_out.write(read)
        else:
            bamp_out.write(read)

    bam_in.close()
    bamp_out.close()
    bamm_out.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
