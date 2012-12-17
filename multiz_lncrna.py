#!/usr/bin/env python
from optparse import OptionParser
from pygr import worldbase
import gff

################################################################################
# multiz_lncrna.py
#
# Return hg19 46-way multiz alignments of a specified lncRNA gene. By default,
# just do the exons, but as an option do the entire span.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gene id>'
    parser = OptionParser(usage)
    parser.add_option('-l', dest='lncrna_gtf', default='/Users/dk/research/common/data/lncrna/lnc_catalog.gtf', help='lncRNA gtf file [Default: %default]')
    parser.add_option('-s', dest='span', action='store_true', default=False, help='Map the gene\'s entire span, i.e. introns too [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide gene id')
    else:
        gene_id = args[0]

    # get human genome
    hg19 = worldbase.Bio.Seq.Genome.HUMAN.hg19()

    # get gene exon intervals
    gene_ivals = []
    for line in open(options.lncrna_gtf):
        a = line.split('\t')
        if gff.gtf_kv(a[8])['gene_id'] == gene_id:
            chrom = a[0]
            start = int(a[3])
            end = int(a[4])
            # ignoring orientation at the moment

            gene_ivals.append(hg19[chrom][start:end])

    # get hg19 msa
    msa = worldbase.Bio.MSA.UCSC.hg19_multiz46way()

    # map returned sequences back to genome name
    idDict = ~(msa.seqDict)

    # print alignments
    for gi in gene_ivals:
        for src, dest, edg in msa[gi].edges():
            print repr(gi), repr(src), repr(dest), idDict[dest], edg.length()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
