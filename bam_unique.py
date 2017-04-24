#!/usr/bin/env python
from optparse import OptionParser
import pysam

################################################################################
# bam_unique.py
#
# Remove multi-mapping alignments from a BAM file
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <input_bam> <output_bam>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input and output BAMs')
    else:
        input_bam = args[0]
        output_bam = args[1]

    bam_in = pysam.Samfile(input_bam, 'rb')
    bam_out = pysam.Samfile(output_bam, 'wb', template=bam_in)

    for aligned_read in bam_in:
        if aligned_read.get_tag('NH') == 1:
            bam_out.write(aligned_read)

    bam_in.close()
    bam_out.close()

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
