#!/usr/bin/env python
from optparse import OptionParser
import pysam

################################################################################
# set_bam_xs.py
#
# Set the XS tag properly in a BAM file. Currently assumes first-strand because
# that's all that I've seen so far.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam>'
    parser = OptionParser(usage)
    parser.add_option('-o', dest='bam_out_file', help='Output BAM file')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide BAM file')
    else:
        bam_in_file = args[0]

    if options.bam_out_file == None:
        b, e = os.path.splitext(bam_in_file)
        options.bam_out_file = b+'_xs'+e

    bam_in = pysam.Samfile(bam_in_file, 'rb')
    bam_out = pysam.Samfile(options.bam_out_file, 'wb', template=bam_in)

    for aligned_read in bam_in:
        # remove existing XS tag
        rm_xs(aligned_read)

        # set XS tag
        if aligned_read.is_paired:
            if aligned_read.is_read1:
                if aligned_read.is_reverse:
                    aligned_read.tags = aligned_read.tags + [('XS','+')]
                else:
                    aligned_read.tags = aligned_read.tags + [('XS','-')]
            else:
                if aligned_read.is_reverse:
                    aligned_read.tags = aligned_read.tags + [('XS','-')]
                else:
                    aligned_read.tags = aligned_read.tags + [('XS','+')]
        else:
            if aligned_read.is_reverse:
                aligned_read.tags = aligned_read.tags + [('XS','+')]
            else:
                aligned_read.tags = aligned_read.tags + [('XS','-')]

        # output
        bam_out.write(aligned_read)
    

################################################################################
# rm_xs
#
# Remove the XS tag from the AlignedRead object
################################################################################
def rm_xs(aligned_read):
    xs_i = 0
    while xs_i < len(aligned_read.tags) and aligned_read.tags[xs_i][0] != 'XS':
        xs_i += 1

    if xs_i < len(aligned_read.tags):
        aligned_read.tags = aligned_read.tags[:xs_i] + aligned_read.tags[xs_i+1:]


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
