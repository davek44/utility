#!/usr/bin/env python
from optparse import OptionParser
import math
import cufflinks, gff, stats

################################################################################
# most_expr_isoform.py
#
# Filter a GTF file to include only the most expressed isoform of each gene
# locus.
#
# If every isoform failed in quantification os is unexpressed, we won't return
# any.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf> <fpkm tracking | diff>'
    parser = OptionParser(usage)
    parser.add_option('-a', dest='all_isoforms', default=False, action='store_true', help='Consider all isoforms. Default is to ignore bs ones')
    parser.add_option('-p', dest='pseudocount', default=0.125)
    parser.add_option('-r', dest='random_zeros', default=False, action='store_true', help='Randomly choose an isoform for zero FPKM genes [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide gtf file and fpkm tracking')
    else:
        gtf_file = args[0]
        fpkm_file = args[1]

    gene_max_iso = map_genes(gtf_file, fpkm_file, options.pseudocount, options.all_isoforms, options.random_zeros)

    # filter gtf file
    for line in open(gtf_file):
        a = line.split('\t')
        kv = gff.gtf_kv(a[8])
        gene_id = kv['gene_id']
        tid = kv['transcript_id']        
        
        if gene_max_iso.get(gene_id,None) == tid:
            print line,


################################################################################
# cuff_fpkm
#
# Hash the mean log2 FPKM of all genes in a .fpkm_tracking file.
################################################################################
def cuff_fpkm(fpkm_file, pseudocount):
    cuff = cufflinks.fpkm_tracking(fpkm_tracking_file)

    gene_fpkm = {}
    for gene_id in cuff.genes:
        gene_fpkm[gene_id] = stats.mean([math.log(pseudocount+e,2) for e in cuff.gene_expr(gene_id, not_found=0, fail=0)])

    return gene_fpkm


################################################################################
# diff_fpkm
#
# Hash the mean log2 FPKM of all genes in a cuffdiff output file.
#
# The code is written for the case where there are only two samples compared,
# but it should be OK even more more.
################################################################################
def diff_fpkm(diff_file, pseudocount):
    gene_fpkms = {}

    diff_in = open(diff_file)
    diff_in.readline()
    for line in diff_in:
        a = line.split('\t')

        gene_id = a[0]
        sample1 = a[4]
        sample2 = a[5]
        status = a[6]
        fpkm1 = float(a[7])
        fpkm2 = float(a[8])

        if status == 'OK':
            if gene_id in gene_fpkms:
                gene_fpkms[gene_id] += [fpkm1, fpkm2]
            else:
                gene_fpkms[gene_id] = [fpkm1, fpkm2]

    diff_in.close()

    gene_fpkm = {}
    for gene_id in gene_fpkms:
        log_fpkms = [math.log(fpkm+pseudocount,2) for fpkm in gene_fpkms[gene_id]]
        gene_fpkm[gene_id] = stats.mean(log_fpkms)

    return gene_fpkm


################################################################################
# map_genes
#
# Input
#  gtf_file:      GTF file describing genes.
#  fpkm_file:     Cufflinks .fpkm_tracking or .diff file for the genes.
#  pseudocount:   Pseudocount used for log2 FPKM.
#  all_isoforms:  Consider all isoforms, or filter crap.
# 
# Output
#  gene_max_iso:  Dict mapping gene_id to transcript_id or None
################################################################################
def map_genes(gtf_file, fpkm_file, pseudocount=0.125, all_isoforms=False, random_zeros=False):
    # get expression data
    if fpkm_file[-5:] == '.diff':
        transcript_fpkm = diff_fpkm(fpkm_file, pseudocount)
    else:
        transcript_fpkm = cuff_fpkm(fpkm_file, pseudocount)

    # get genes
    if all_isoforms:
        g2t = gff.g2t(gtf_file)
    else:
        g2t = {}
        for line in open(gtf_file):
            a = line.split('\t')
            kv = gff.gtf_kv(a[8])

            if kv['transcript_type'] not in ['intron', 'prerna', 'nonsense_mediated_decay', 'retained_intron', 'non_stop_decay']:
                g2t.setdefault(kv['gene_id'],set()).add(kv['transcript_id'])

    # map gene_id's to max expression isoform
    gene_max_iso = {}
    min_fpkm = math.log(pseudocount, 2)
    for gid in g2t:
        max_fpkm_tid = None
        max_fpkm = min_fpkm

        for tid in g2t[gid]:
            if transcript_fpkm.get(tid,min_fpkm) > max_fpkm:
                max_fpkm_tid = tid
                max_fpkm = transcript_fpkm[tid]

        gene_max_iso[gid] = max_fpkm_tid

    # choose isoforms for None
    if random_zeros:
        for gid in g2t:
            if gene_max_iso[gid] == None:
                gene_max_iso[gid] = random.choice(g2t[gid])

    return gene_max_iso


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
