#!/usr/bin/env python
from optparse import OptionParser
import pysam

'''
bam_quality.py

Remove low quality alignments from a BAM file.
'''


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <input_bam> <output_bam>'
    parser = OptionParser(usage)
    parser.add_option('-q', dest='mapq_t',
        type='int', default=None,
        help='Keep alignments with mapping quality at or above [Default: %default]')
    parser.add_option('-s', dest='score_t',
        type='int', default=None,
        help='Keep alignments with alignment score at or above [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input and output BAMs')
    else:
        input_bam = args[0]
        output_bam = args[1]

    bam_in = pysam.AlignmentFile(input_bam, 'r')
    bam_out = pysam.AlignmentFile(output_bam, 'wb', template=bam_in)

    for align in bam_in:
        if options.mapq_t is None or align.mapping_quality >= options.mapq_t:
            if options.score_t is None or align.get_tag('AS') >= options.score_t:
                bam_out.write(align)

    bam_in.close()
    bam_out.close()

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
