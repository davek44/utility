#!/usr/bin/env python
from optparse import OptionParser
import gff, util
import os, subprocess, sys

################################################################################
# gtf_multimaps.py
#
# Print a summary table about multimapping reads for the transcripts in a gtf
# file.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf file> <bam file>'
    parser = OptionParser(usage)
    parser.add_option('-i', dest='intersect_done', default=False, action='store_true', help='intersectBed is already done [Default: %default]')
    parser.add_option('-o', dest='output_prefix', help='Prefix for the intersectBed intermediate file [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide gtf file and bam file')
    else:
        gtf_file = args[0]
        bam_file = args[1]

    if options.output_prefix:
        ib_file = '%s_reads_genes.gff' % options.output_prefix
    else:
        ib_file = 'reads_genes.gff'

    if not options.intersect_done:
        # overlap genes w/ aligned reads
        p = subprocess.Popen('intersectBed -s -wo -abam -bed -a %s -b %s > %s' % (bam_file,gtf_file,ib_file), shell=True)
        os.waitpid(p.pid,0)

    # count transcriptome alignments per read
    read_aligns = {}
    for line in open(ib_file):
        a = line.split('\t')
        chrom = a[0]
        start = int(a[1])
        read_id = a[3]

        read_aligns.setdefault(read_id,set()).add((chrom,start))

    # hash reads by gene
    gene_reads = {}
    for line in open(ib_file):
        a = line.split('\t')
        read_id = a[3]
        gene_id = gff.gtf_kv(a[14])['transcript_id']
        gene_reads.setdefault(gene_id,[]).append(read_id)

    # print gene stats
    for gene_id in gene_reads:
        align_counts = [len(read_aligns[read_id]) for read_id in gene_reads[gene_id]]
        multi_count = float(len([ac for ac in align_counts if ac > 1]))
        cols = (gene_id, len(align_counts), util.mean(align_counts), multi_count/float(len(align_counts)))
        print '%-15s %7d %7.2f %7.2f' % cols


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
