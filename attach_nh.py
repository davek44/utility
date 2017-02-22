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
    write_header(sam_in.header)

    sam_out = pysam.AlignmentFile('-', 'w', template=sam_in)

    last_id = 'not a header'
    read1_aligns = []
    read2_aligns = []

    for align in sam_in:
        if align.is_unmapped:
            # read stream concludes
            output_read(sam_out, read1_aligns)
            output_read(sam_out, read2_aligns)
            read1_aligns = []
            read2_aligns = []

        else:
            read_id = align.query_name

            if not match_id(read_id, last_id, align.is_paired):
                # read stream concludes

                # output
                output_read(sam_out, read1_aligns)
                output_read(sam_out, read2_aligns)

                # reset
                read1_aligns = []
                read2_aligns = []

            # read stream continues
            if align.is_read1:
                read1_aligns.append(align)
            else:
                read2_aligns.append(align)

            # update read id
            last_id = read_id

    sam_in.close()
    sam_out.close()


def match_id(id1, id2, paired):
    ''' Match read_id's.

    First case handles most datasets.
    Second case handles paired end datasets where they got fancy.
    '''

    return id1 == id2 or (paired and id1[:-1] == id2[:-1] and id1[-1] in '12' and id2[-1] in '12')


def output_read(sam_out, read_aligns):
    nh_tag = len(read_aligns)
    for align in read_aligns:
        align.set_tag('NH',nh_tag)
        sam_out.write(align)


def write_header(header):
    hd = header['HD']
    print('@HD\tVN:%s\tSO:%s' % (hd['VN'], hd['SO']))

    for sq in header['SQ']:
        print('@SQ\tSN:%s\tLN:%d' % (sq['SN'],sq['LN']))

    pg = header['PG'][0]
    print('@PG\tID:%s\tPN:%s\tVN:%s\tCL:%s' % (pg['ID'],pg['PN'],pg['VN'],pg['CL']), flush=True)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
