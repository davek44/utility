#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess
import networkx
import gff

################################################################################
# clean_gtf.py
#
# Clean up a gtf file (like the messy RefSeq download) to match the linc catalog.
#
# That consists of:
# 1. Sort out duplicated(?) genes that have the same id's and make them unique
# 2. Try to assign isoforms from the same gene a similar gene id
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error(usage)
    else:
        gtf_file = args[0]

    ############################################
    # fix multi-chromosome genes
    ############################################
    # find multi-chromosome genes
    tx_chrs = {}
    for line in open(gtf_file):
        a = line.split('\t')
        kv = gff.gtf_kv(a[8])
        if kv['transcript_id'].startswith('NM_'):
            tx_chrs.setdefault(kv['transcript_id'],set()).add(a[0])
    multi_genes = set([tid for tid in tx_chrs if len(tx_chrs[tid]) > 1])

    # revise gtf
    tx_gene = {}
    gtf_out = open('tmp.gtf','w')
    for line in open(gtf_file):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()
        kv = gff.gtf_kv(a[8])

        # if multi-chrom gene, supplement id's
        if kv['transcript_id'] in multi_genes:
            kv['transcript_id'] += 'c%s'%a[0][3:]
            a[8] = gff.kv_gtf(kv)

        # map trans to gene (forget the actual gene id's; they don't consider "_dup")
        tx_gene[kv['transcript_id']] = kv['transcript_id']

        # print new line
        print >> gtf_out, '\t'.join(a)

    gtf_out.close()

    ############################################
    # merge transcripts into genes
    ############################################
    # intersect and build overlapping transcript graph
    G = networkx.Graph()
    p = subprocess.Popen('intersectBed -f 0.2 -r -wo -s -a tmp.gtf -b tmp.gtf', shell=True, stdout=subprocess.PIPE)
    line = p.stdout.readline()
    while line:
        a = line.split('\t')
        tid1 = gff.gtf_kv(a[8])['transcript_id']
        tid2 = gff.gtf_kv(a[17])['transcript_id']
        G.add_edge(tid1,tid2)
        line = p.stdout.readline()
    p.communicate()

    # combine connected components as genes
    for component in networkx.algorithms.components.connected.connected_components(G):
        comp_gene = 'G'+tx_gene[component[0]]        
        for tid in component:
            tx_gene[tid] = comp_gene
    for tid in tx_gene:
        if tx_gene[tid][0] != 'G':
            tx_gene[tid] = 'G'+tx_gene[tid]

    ############################################
    # output
    ############################################
    # print
    for line in open('tmp.gtf'):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        kv = gff.gtf_kv(a[8])
        kv['gene_id'] = tx_gene[kv['transcript_id']]
        a[8] = gff.kv_gtf(kv)
        
        print '\t'.join(a)

    # clean
    os.remove('tmp.gtf')


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
