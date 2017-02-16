#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
import math, sys
import cufflinks, gff, stats

################################################################################
# promoters.py
#
# Return promoter regions in GFF format from a reference GTF.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <ref_gtf>'
    parser = OptionParser(usage)
    #parser.add_option()
    parser.add_option('-d', dest='downstream', type='int', default=1000, help='Downstream bp for promoters [Default: %default]')
    parser.add_option('-f', dest='fpkm_tracking', help='Use cufflinks FPKM estimates to choose the most expressed isoform')
    parser.add_option('-u', dest='upstream', type='int', default=1000, help='Upstream bp for promoters [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide reference GTF')
    else:
        ref_gtf = args[0]

    g2t = gff.g2t(ref_gtf)
    transcripts = gff.read_genes(ref_gtf)
    source = open(ref_gtf).readline().split()[1]

    if options.fpkm_tracking:
        iso_fpkm_tracking = cufflinks.fpkm_tracking(options.fpkm_tracking)

    for gene_id in g2t:
        gene_transcripts = list(g2t[gene_id])
        gene_strand = transcripts[gene_transcripts[0]].strand
        if gene_strand not in ['+','-']:
            print('WARNING: %s discluded for lack of strand' % gene_id, file=sys.stderr)
            continue

        # choose TSS
        if options.fpkm_tracking:
            # find most expressed isoform
            promoter_tid = gene_transcripts[0]
            max_fpkm = stats.geo_mean([1+fpkm for fpkm in iso_fpkm_tracking.gene_expr(promoter_tid)])
            for transcript_id in gene_transcripts[1:]:
                transcript_fpkm = stats.geo_mean([1+fpkm for fpkm in iso_fpkm_tracking.gene_expr(transcript_id)])
                if math.isnan(max_fpkm) or transcript_fpkm > max_fpkm:
                    promoter_tid = transcript_id
                    max_fpkm = transcript_fpkm

            # get isoform tss
            if gene_strand == '+':
                tss = transcripts[promoter_tid].exons[0].start
            else:
                tss = transcripts[promoter_tid].exons[-1].end

        else:
            # find most upstream tss
            promoter_tid = gene_transcripts[0]
            if gene_strand == '+':
                upstream_tss = transcripts[promoter_tid].exons[0].start
            else:
                upstream_tss = transcripts[promoter_tid].exons[-1].end

            for transcript_id in gene_transcripts[1:]:
                if gene_strand == '+':
                    transcript_pos = transcripts[transcript_id].exons[0].start
                    if transcript_pos < upstream_tss:
                        promoter_tid = transcript_id
                        upstream_tss = transcript_pos
                else:
                    transcript_pos = transcripts[transcript_id].exons[-1].end
                    if transcript_pos > upstream_tss:
                        promoter_tid = transcript_id
                        upstream_tss = transcript_pos

            tss = upstream_tss

        # print promoter from the tss
        if gene_strand == '+':
            if tss - options.upstream < 1:
                print('WARNING: %s discluded for nearness to chromosome end' % gene_id, file=sys.stderr)
            else:
                tx = transcripts[promoter_tid]
                cols = [tx.chrom, source, 'promoter', str(tss-options.upstream), str(tss+options.downstream), '.', tx.strand, '.', gff.kv_gtf(tx.kv)]
                print('\t'.join(cols))

        else:
            if tss - options.downstream < 1:
                print('WARNING: %s discluded for nearness to chromosome end' % gene_id, file=sys.stderr)
            else:
                tx = transcripts[promoter_tid]
                cols = [tx.chrom, source, 'promoter', str(tss-options.downstream), str(tss+options.upstream), '.', tx.strand, '.', gff.kv_gtf(tx.kv)]
                print('\t'.join(cols))
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
