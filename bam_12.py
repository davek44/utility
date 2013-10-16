#!/usr/bin/env python
from optparse import OptionParser
import pdb, os
import pysam

################################################################################
# bam_12.py
#
# Separate the alignments in a BAM file into two BAM files of the first and
# second reads.
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

    bam_pre = os.path.splitext(bam_file)[0]

    bam_in = pysam.Samfile(bam_file, 'rb')
    bam1_out = pysam.Samfile('%s_1.bam'%bam_pre, 'wb', header=bam_in.header)
    bam2_out = pysam.Samfile('%s_2.bam'%bam_pre, 'wb', header=bam_in.header)

    for read in bam_in:
        if read.is_read1:
            bam1_out.write(read)
        else:
            bam2_out.write(read)

    bam_in.close()
    bam1_out.close()
    bam2_out.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
