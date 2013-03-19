#!/usr/bin/env python
from optparse import OptionParser
import math
import cufflinks, gff

################################################################################
# most_expr_isoform.py
#
# Filter a GTF file to include only the most expressed isoform of each gene
# locus.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf> <fpkm tracking>'
    parser = OptionParser(usage)
    parser.add_option('-p', dest='pseudocount', default=0.1)
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide gtf file and fpkm tracking')
    else:
        gtf_file = args[0]
        fpkm_tracking_file = args[1]

    # get expression data
    cuff = cufflinks.fpkm_tracking(fpkm_tracking_file)

    # get genes
    g2t = gff.g2t(gtf_file)

    # map gene_id's to max expression isoform
    gene_max_iso = {}
    for gid in g2t:
        max_expr_tid = None
        for tid in g2t[gid]:
            expr_mean = sum([math.log(options.pseudocount+e,2) for e in cuff.gene_expr(tid)])
            if max_expr_tid == None or expr_mean > max_expr:
                max_expr = expr_mean
                max_expr_tid = tid

        gene_max_iso[gid] = max_expr_tid

    # filter gtf file
    for line in open(gtf_file):
        a = line.split('\t')
        kv = gff.gtf_kv(a[8])
        gene_id = kv['gene_id']
        tid = kv['transcript_id']        
        
        if gene_max_iso[gene_id] == tid:
            print line,
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
