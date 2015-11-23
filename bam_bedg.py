#!/usr/bin/env python
from optparse import OptionParser
import pysam

################################################################################
# bam_bedg
#
# Map a BAM file of aligned reads from a ChIP-seq or ATAC-seq to a BEDGRAPH
# file, counting only the events relevant to that experiment.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam> <bedg>'
    parser = OptionParser(usage)
    parser.add_option('-l', dest='frag_len', default=200, type='int', action='store_true')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input BAM and output BEDGRAPH files')
    else:
        bam_file = args[0]
        bedg_file = args[1]

    chrom_events = {}

    bam_in = pysam.Samfile(bam_file, 'rb')
    for align in bam_in:
        # get chrom
        chrom = bam_in.references[align.tid]

        # weight multi-mappers
        multi_weight = weight_multi(align)

        # determine fragment length
        if align.is_proper_pair:
            frag_len = abs(align.tlen)
        else:
            frag_len = options.frag_len

        # map to event position
        event_pos = align.reference_start + frag_len/2

        # save
        if chrom not in chrom_events:
            chrom_events[chrom] = {}
        chrom_events[chrom][event_pos] = chrom_events[chrom].get(event_pos,0) + multi_weight
    bam_in.close()

    # output BEDGRAPH


def weight_multi(align):
    ''' Weight the alignment by its multimap properties

    I'm making this a separate function, because I might
    want to use more sophisticated weights later.
    '''
    try:
        nh_tag = align_read.opt('NH')
    except:
        nh_tag = 1

    multi_weight = 1.0 / nh_tag

    return multi_weight



################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
