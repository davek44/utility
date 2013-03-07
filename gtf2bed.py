#!/usr/bin/env python
from optparse import OptionParser
import gff

################################################################################
# gtf2bed.py
#
# Convert a gtf file to a bed file.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf file>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='cds', action='store_true', default=False, help='Use CDS, not exons [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide gtf file')
    else:
        gtf_file = args[0]

    genes = gff.read_genes(gtf_file)

    for transcript_id in genes:
        g = genes[transcript_id]

        if options.cds:
            block_sizes = ','.join([str(ex.end-ex.start+1) for ex in g.cds])
            block_starts = ','.join([str(ex.start-g.cds[0].start) for ex in g.cds])

            cols = [g.chrom, str(g.cds[0].start-1), str(g.cds[-1].end), transcript_id, '0', g.strand, '0', '0', '255,0,0', str(len(g.cds)), block_sizes, block_starts]

        else:
            block_sizes = ','.join([str(ex.end-ex.start+1) for ex in g.exons])
            block_starts = ','.join([str(ex.start-g.exons[0].start) for ex in g.exons])

            cols = [g.chrom, str(g.exons[0].start-1), str(g.exons[-1].end), transcript_id, '0', g.strand, '0', '0', '255,0,0', str(len(g.exons)), block_sizes, block_starts]

        print '\t'.join(cols)
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
