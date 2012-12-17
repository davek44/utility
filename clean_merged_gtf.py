#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess, sys
import networkx as nx
import gff

################################################################################
# clean_merged_gtf.py
#
# Disconnect reference genes that were combined by a new overlapping transcript
# and make sure all loci overlapping a reference gene are named by that
# gene id and it's short gene name.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <ref gtf> <merged gtf>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error(usage)
    else:
        ref_gtf = args[0]
        merged_gtf = args[1]

    # get mappings
    ref_t2g = gff.t2g(ref_gtf)
    merged_t2g = gff.t2g(merged_gtf)
    merged_g2t = gff.g2t(merged_gtf)

    # hash gene_name's by tid
    ref_gid_names = {}
    for line in open(ref_gtf):
        a = line.split('\t')
        kv = gff.gtf_kv(a[8])
        if 'gene_name' in kv:
            ref_gid_names[kv['gene_id']] = kv['gene_name']

    # hash merged lines by tid
    merged_tid_lines = {}
    for line in open(merged_gtf):
        a = line.split('\t')
        tid = gff.gtf_kv(a[8])['transcript_id']
        merged_tid_lines.setdefault(tid,[]).append(line)

    # intialize orphan gene_id
    orphan_num = 1

    for mgene_id in merged_g2t:
        # count reference genes
        ref_genes = set()
        for tid in merged_g2t[mgene_id]:
            if tid in ref_t2g:
                ref_genes.add(ref_t2g[tid])

        # if no known genes, leave it alone
        if len(ref_genes) == 0:
            for tid in merged_g2t[mgene_id]:
                print ''.join(merged_tid_lines[tid]),

        # if known gene, set gene_id to it
        elif len(ref_genes) == 1:
            new_gene_id = list(ref_genes)[0]
            for tid in merged_g2t[mgene_id]:
                for line in merged_tid_lines[tid]:
                    a = line.split('\t')
                    kv = gff.gtf_kv(a[8])
                    kv['gene_id'] = new_gene_id
                    if new_gene_id in ref_gid_names:
                        kv['gene_name'] = ref_gid_names[new_gene_id]
                    a[8] = gff.kv_gtf(kv)
                    print '\t'.join(a)

        # if two known genes were combined, fix it
        elif len(ref_genes) > 1:
            # compute transcript overlaps and build overlap graph
            tid_overlap_graph = make_overlap_graph(mgene_id, merged_g2t, merged_tid_lines)

            # map each new transcript to the ref gene_id's overlapped
            tid_ref_genes = {}
            for (tid1,tid2) in tid_overlap_graph.edges():
                if tid1 in ref_t2g and tid2 not in ref_t2g:
                    tid_ref_genes.setdefault(tid2,set()).add(ref_t2g[tid1])
                elif tid1 not in ref_t2g and tid2 in ref_t2g:
                    tid_ref_genes.setdefault(tid1,set()).add(ref_t2g[tid2])

            # remove new transcripts overlapping multiple ref gene_id's
            for tid in tid_ref_genes:
                if len(tid_ref_genes[tid]) > 1:
                    print >> sys.stderr, 'Removing %s' % tid
                    tid_overlap_graph.remove_node(tid)

            # remove edges connecting separate reference genes
            for (tid1,tid2) in tid_overlap_graph.edges():
                if tid1 in ref_t2g and tid2 in ref_t2g and ref_t2g[tid1] != ref_t2g[tid2]:
                    tid_overlap_graph.remove_edge(tid1,tid2)

            # map to new gene_id's; missing means eliminate transcript
            tid_new_gid, orphan_num = map_new_gid(tid_overlap_graph, orphan_num, ref_t2g)

            for tid in merged_g2t[mgene_id]:
                if tid in tid_new_gid:
                    for line in merged_tid_lines[tid]:
                        a = line.split('\t')
                        kv = gff.gtf_kv(a[8])
                        kv['gene_id'] = tid_new_gid[tid]
                        if tid_new_gid[tid] in ref_gid_names:
                            kv['gene_name'] = ref_gid_names[tid_new_gid[tid]]
                        a[8] = gff.kv_gtf(kv)
                        print '\t'.join(a)


################################################################################
# make_overlap_graph
#
# Compute transcript overlaps and build overlap graph
################################################################################
def make_overlap_graph(mgene_id, merged_g2t, merged_tid_lines):
    # make temporary gff file for gene
    tmp_out = open('%s.gff' % mgene_id, 'w')
    for tid in merged_g2t[mgene_id]:
        for line in merged_tid_lines[tid]:
            a = line.split('\t')
            if a[2] == 'exon':
                print >> tmp_out, line,
    tmp_out.close()

    tid_overlap_graph = nx.Graph()

    # intersect with self
    proc = subprocess.Popen('intersectBed -wo -s -a %s.gff -b %s.gff' % (mgene_id,mgene_id), shell=True, stdout=subprocess.PIPE)
    line = proc.stdout.readline()
    while line:
        a = line.split('\t')

        tid1 = gff.gtf_kv(a[8])['transcript_id']
        tid2 = gff.gtf_kv(a[17])['transcript_id']

        # ignore same and ignore different ref genes
        if tid1 != tid2:                    
            tid_overlap_graph.add_edge(tid1,tid2)

        line = proc.stdout.readline()
    proc.communicate()

    os.remove('%s.gff' % mgene_id)

    return tid_overlap_graph


################################################################################
# map_new_gid
#
# Map to new gene_id's; None means eliminate the transcript
################################################################################
def map_new_gid(tid_overlap_graph, orphan_num, ref_t2g):
    tid_new_gid = {}

    # for each connected component
    for cc in nx.connected_component_subgraphs(tid_overlap_graph):
        # assign the reference gene id
        cc_ref_genes = set()
        for tid in cc.nodes():
            if tid in ref_t2g:
                cc_ref_genes.add(ref_t2g[tid])

        if len(cc_ref_genes) == 0:
            # orphaned: assign new gene_id
            new_gene_id = 'orphan_XLOC_%06d' % orphan_num
            orphan_num += 1
            for tid in cc.nodes():
                tid_new_gid[tid] = new_gene_id

        elif len(cc_ref_genes) == 1:
            # set all to that gene_id
            new_gene_id = list(cc_ref_genes)[0]
            for tid in cc.nodes():
                tid_new_gid[tid] = new_gene_id

        else:                    
            print >> sys.stderr, 'Ambiguous connected component: %s' % (' '.join(list(cc_ref_genes)))
            # blow away all new transcripts
            # set references ones to their old gene_id
            for tid in cc.nodes():
                if tid in ref_t2g:
                    tid_new_gid[tid] = ref_t2g[tid]

    return tid_new_gid, orphan_num
                    
################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
