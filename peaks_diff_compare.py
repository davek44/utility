#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess
import ggplot, gff, math, os, stats, subprocess, tempfile
import ripseq

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

from scipy.stats import hypergeom

################################################################################
# peaks_diff_compare.py
#
# Compare RNAs bound according to some peak calls to a cuffdiff run.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <peaks gff> <diff>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='clip_fpkm_file', help='Control FPKM tracking file')
    parser.add_option('-g', dest='ref_gtf', default='%s/gencode.v18.annotation.gtf'%os.environ['GENCODE'])
    parser.add_option('--ggplot', dest='ggplot_script', default='%s/peaks_diff_compare.r'%os.environ['RDIR'], help='Script to make plots with [Default: %default]')
    parser.add_option('-m', dest='max_stat', default=10, type='float', help='Max cuffdiff stat [Default: %default]')
    parser.add_option('-o', dest='output_pre', default='', help='Output prefix [Default: %default]')
    parser.add_option('-r', dest='rbp', default='RBP', help='RBP name [Default: %default]')
    parser.add_option('-s', dest='single_gene_loci', default=False, action='store_true', help='Only use single gene loci [Default: %default]')
    parser.add_option('-t', dest='test_stat', default=False, action='store_true', help='Use test statistic rather than fold change [Default: %default]')
    parser.add_option('--sample1', dest='sample1', help='Sample_1 name in cuffdiff')
    parser.add_option('--sample2', dest='sample2', help='Sample_2 name in cuffdiff')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide peaks GFF and .diff file')
    else:
        peaks_gff = args[0]
        diff_file = args[1]

    ##################################################
    # process GTF
    ##################################################
    if options.single_gene_loci:
        single_gtf_fd, single_gtf_file = filter_single(options.ref_gtf)
        options.ref_gtf = single_gtf_file

    gtf_genes = gff.gtf_gene_set(options.ref_gtf)

    ##################################################
    # collect CLIP peak bound genes
    ##################################################
    peak_genes = set()
    p = subprocess.Popen('intersectBed -s -u -a %s -b %s' % (options.ref_gtf, peaks_gff), shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        peak_genes.add(gff.gtf_kv(line.split('\t')[8])['gene_id'])
    p.communicate()

    # find expressed genes in peak calls
    silent_genes = set()
    if options.clip_fpkm_file:
        silent_genes = find_silent(options.clip_fpkm_file)

    ##################################################
    # collect RIP stats
    ##################################################
    if options.test_stat:
        rip_fold, rip_bound = ripseq.hash_rip(diff_file, just_ok=True, use_fold=False, max_stat=options.max_stat, one_rbp=True)
    else:
        rip_fold, rip_bound = ripseq.hash_rip(diff_file, use_fold=True, max_stat=options.max_stat, one_rbp=True)
        rip_fold = ripseq.hash_rip_fold(diff_file, min_fpkm=0.125, max_fold=10, one_rbp=True)

    # TEMP: print bound genes
    # genes_out = open('%s_genes.txt' % options.output_pre, 'w')
    # for gene_id in rip_bound:
    #     if rip_bound[gene_id]:
    #         print >> genes_out, gene_id, rip_fold[gene_id]
    # genes_out.close()

    ##################################################
    # plot bound and unbound distributions
    ##################################################
    # construct data frame
    df_dict = {'Gene':[], 'CLIP':[], 'RIP':[]}
    for gene_id in rip_fold:
        if gene_id in gtf_genes and (len(silent_genes) == 0 or gene_id not in silent_genes):
            df_dict['Gene'].append(gene_id)
            df_dict['RIP'].append(rip_fold[gene_id])
            if gene_id in peak_genes:
                df_dict['CLIP'].append('Bound')
            else:
                df_dict['CLIP'].append('Unbound')

    ggplot.plot(options.ggplot_script, df_dict, [options.output_pre, options.rbp, options.test_stat])

    ##################################################
    # compute stats on bound and unbound distributions
    ##################################################
    bound_fold = [df_dict['RIP'][i] for i in range(len(df_dict['RIP'])) if df_dict['CLIP'][i] == 'Bound']
    unbound_fold = [df_dict['RIP'][i] for i in range(len(df_dict['RIP'])) if df_dict['CLIP'][i] == 'Unbound']

    # perform statistical test
    z, p = stats.mannwhitneyu(bound_fold, unbound_fold)

    stats_out = open('%s_stats.txt' % options.output_pre, 'w')
    cols = (options.rbp, len(bound_fold), stats.mean(bound_fold), len(unbound_fold), stats.mean(unbound_fold), z, p)
    print >> stats_out, '%-10s  %5d  %6.2f  %5d  %6.2f  %6.2f  %9.2e' % cols
    stats_out.close()

    ##################################################
    # plot venn diagram
    ##################################################
    rip_genes = set([df_dict['Gene'][i] for i in range(len(df_dict['Gene'])) if rip_bound.get(df_dict['Gene'][i],False)])

    clip_only = len(peak_genes - rip_genes)
    rip_only = len(rip_genes - peak_genes)
    both = len(peak_genes & rip_genes)

    if options.clip_fpkm_file:
        print >> sys.stderr, 'Ignoring silent genes for hypergeometric test'

    # k is x
    # K is n
    # N is M
    # n is N
    # hypergeom.sf(x, M, n, N, loc=0)

    p1 = hypergeom.sf(both-1, len(gtf_genes), len(peak_genes), len(rip_genes))
    p2 = hypergeom.sf(both-1, len(gtf_genes), len(rip_genes), len(peak_genes))

    hyper_out = open('%s_hyper.txt' % options.output_pre, 'w')
    cols = (p1, p2, both, clip_only, rip_only, len(peak_genes), len(rip_genes), len(gtf_genes))
    print >> hyper_out, '%7.2e  %7.2e  %5d  %5d  %5d  %5d  %5d %5d' % cols
    hyper_out.close()

    if clip_only > 0 and rip_only > 0:
        plt.figure()
        # venn_diag = venn2(subsets=(clip_only, rip_only, both), set_labels=['CLIP', 'fRIP'], set_colors=['#e41a1c', '#377eb8'])
        # venn_diag = venn2(subsets=(clip_only, rip_only, both), set_labels=['CLIP', 'fRIP'], set_colors=['#e41a1c', '#1ae47d'])
        venn_diag = venn2(subsets=(clip_only, rip_only, both), set_labels=['CLIP', 'fRIP'], set_colors=['#e41a1c', '#A1A838'])
        plt.savefig('%s_venn.pdf' % options.output_pre)

    ##################################################
    # clean
    ##################################################
    if options.single_gene_loci:
        os.close(single_gtf_fd)
        os.remove(single_gtf_file)


################################################################################
# filter_single
#
# Input
#   ref_gtf:
#
# Output
#   single_gtf_fd:
#   single_gtf_file:
################################################################################
def filter_single(ref_gtf):
    # intersect with self and compute overlap sets
    #p = subprocess.Popen('intersectBed -sorted -wo -s -a %s -b %s' % (ref_gtf, ref_gtf), shell=True, stdout=subprocess.PIPE)
    p = subprocess.Popen('intersectBed -wo -s -a %s -b %s' % (ref_gtf, ref_gtf), shell=True, stdout=subprocess.PIPE)

    # computer overlaps
    gene_overlaps = {}
    for line in p.stdout:
        a = line.split('\t')

        gid1 = gff.gtf_kv(a[8])['gene_id']
        gid2 = gff.gtf_kv(a[17])['gene_id']

        if gid1 != gid2:
            gene_overlaps.setdefault(gid1,set()).add(gid2)
            gene_overlaps.setdefault(gid2,set()).add(gid1)

    p.communicate()

    # filter overlapping genes out
    single_gtf_fd, single_gtf_file = tempfile.mkstemp()
    single_gtf_out = open(single_gtf_file, 'w')
    for line in open(ref_gtf):
        a = line.split('\t')
        gene_id = gff.gtf_kv(a[8])['gene_id']
        if gene_id not in gene_overlaps:
            print >> single_gtf_out, line,
    single_gtf_out.close()

    return single_gtf_fd, single_gtf_file


################################################################################
# find_silent
#
# Input:
#  control_fpkm_file: Cufflinks FPKM file.
#  silent_fpkm:       FPKM threshold to call a gene silent.
#
# Output:
#  silent_genes:      Set of silent gene_id's.
################################################################################
def find_silent(clip_fpkm_file, silent_fpkm=0.1):
    # get fpkms (possibly from an isoform file)
    gene_fpkms = {}
    control_fpkm_in = open(clip_fpkm_file)
    control_fpkm_in.readline()
    for line in control_fpkm_in:
        a = line.split('\t')
        gene_id = a[3]
        fpkm = float(a[9])
        gene_fpkms[gene_id] = gene_fpkms.get(gene_id,0) + fpkm
    control_fpkm_in.close()

    silent_genes = set()
    for gene_id in gene_fpkms:
        if gene_fpkms[gene_id] < silent_fpkm:
            silent_genes.add(gene_id)

    return silent_genes


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
