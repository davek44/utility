#!/usr/bin/env python
from optparse import OptionParser
import sys
import pysam

################################################################################
# attach_nh.py
#
# Attach NH tags to a stream of SAM alignments from Bowtie2.
#
# Note: I'm not sure how paired end reads will stream in.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    sam_in = pysam.AlignmentFile('-', 'r')
    sam_out = pysam.AlignmentFile('-', 'w', template=sam_in)

    last_id = None
    read_aligns = []

    for align in sam_in:
        if align.is_unmapped:
            # read stream concludes
            output_read(sam_out, read_aligns)
            read_aligns = []

        else:
            read_id = (align.query_name, align.is_read1)

            if read_id == last_id:
                # read stream continues
                read_aligns.append(align)
            else:
                # read stream concludes
                output_read(sam_out, read_aligns)
                read_aligns = []

            # update read id
            last_id = read_id


def output_read(sam_out, read_aligns):
    nh_tag = len(read_aligns)
    for align in read_aligns:
        align.set_tag('NH',nh_tag)
        sam_out.write(align)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
