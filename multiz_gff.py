#!/usr/bin/env python
from optparse import OptionParser
from pygr import worldbase
import gff

################################################################################
# multiz_gff.py
#
# Return hg19 46-way multiz alignments of the entries in a gff file.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gff file>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide gff file.')
    else:
        gff_file = args[0]

    # get human genome
    hg19 = worldbase.Bio.Seq.Genome.HUMAN.hg19()

    # get feature intervals
    feat_ivals = []
    for line in open(gff_file):
        a = line.split('\t')

        chrom = a[0]
        start = int(a[3])
        end = int(a[4])
        # ignoring orientation at the moment

        feat_ivals.append(hg19[chrom][start:end])

    # get hg19 msa
    msa = worldbase.Bio.MSA.UCSC.hg19_multiz46way()

    # map returned sequences back to genome name
    idDict = ~(msa.seqDict)

    # print alignments
    for gi in feat_ivals:
        for src, dest, edg in msa[gi].edges():
            print repr(gi), repr(src), repr(dest), idDict[dest], edg.length()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
