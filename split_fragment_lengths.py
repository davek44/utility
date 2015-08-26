#!/usr/bin/env python
from optparse import OptionParser
import pysam

################################################################################
# split_fragment_lengths.py
#
# Split a BAM file based on a fragment length threshold.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam> <split_len> <out_pre>'
    parser = OptionParser(usage)
    parser.add_option('-m', dest='max_length', default=1000, help='Threshold length beyond which we ignore the read [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 3:
        parser.error('Must provide BAM file, split length, and output prefix')
    else:
        bam_file = args[0]
        split_len = int(args[1])
        out_pre = args[2]

    bam_in = pysam.Samfile(bam_file, 'rb')

    minus_out = pysam.Samfile('%s_%d-.bam' % (out_pre,split_len), 'wb', template=bam_in)
    plus_out = pysam.Samfile('%s_%d+.bam' % (out_pre,split_len), 'wb', template=bam_in)

    for alignment in bam_in:
        tl = abs(alignment.template_length)
        if tl == 0:
            pass
        elif tl < split_len:
            minus_out.write(alignment)
        elif tl <= options.max_length:
            plus_out.write(alignment)
        else:
            pass

    minus_out.close()
    plus_out.close()

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
