#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import spearmanr, pearsonr
from collections import Counter
import math, os, subprocess, pdb, shutil, sys
import cufflinks, fdr, gff, stats

from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.lib.ggplot2 as ggplot2
grdevices = importr('grDevices')

################################################################################
# cuff_scatter.py
#
# Print a beautiful scatter plot using a cuffdiff file.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <.diff>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='genes_gtf', help='Print only genes in the given GTF file')
    parser.add_option('-p', dest='pseudocount', type='float', default=0.125, help='FPKM pseudocount for taking logs [Default: %default]')
    parser.add_option('-o', dest='out_dir_pre', default='genes', help='Prefix for output directories [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error(usage)
    else:
        diff_file = args[0]
    
    # get gene_ids
    gene_set = set()
    if options.genes_gtf:
        for line in open(options.genes_gtf):
            a = line.split('\t')
            gid = gff.gtf_kv(a[8])['gene_id']
            gene_set.add(gid)

    # initialize diff data structures
    gene_diffs = {}
    gene_fpkm1 = {}
    gene_fpkm2 = {}
    gene_qvals = {}

    # read diff file
    diff_in = open(diff_file)
    headers = diff_in.readline()
    line = diff_in.readline()
    while line:
        a = line.split('\t')

        gene_id = a[1]
        cond1 = a[4]
        cond2 = a[5]
        status = a[6]
        fpkm1 = float(a[7])
        fpkm2 = float(a[8])
        fold_change = float(a[9])
        test_stat = float(a[10])
        qval = float(a[11])

        #if status == 'OK' and abs(test_stat) < 1e6 and (len(gene_set) == 0 or gene_id in gene_set):
        if abs(test_stat) < 1e6 and (len(gene_set) == 0 or gene_id in gene_set):
            gene_diffs.setdefault((cond1,cond2),[]).append(test_stat)
            gene_qvals.setdefault((cond1,cond2),[]).append(qval)
            gene_fpkm1.setdefault((cond1,cond2),[]).append(fpkm1)
            gene_fpkm2.setdefault((cond1,cond2),[]).append(fpkm2)

        line = diff_in.readline()
    diff_in.close()

    # clean plot directory
    if os.path.isdir('%s_scatter' % options.out_dir_pre):
        shutil.rmtree('%s_scatter' % options.out_dir_pre)
    os.mkdir('%s_scatter' % options.out_dir_pre)

    for cond_key in gene_diffs:
        cond1, cond2 = cond_key

        # title plot
        plot_title = '%s vs %s' % (cond1,cond2)

        # set statistic range        
        #stat_min = stats.quantile(gene_diffs[cond_key], .005)
        #stat_max = stats.quantile(gene_diffs[cond_key], .995)

        # construct data frame
        fpkm1_r = ro.FloatVector([math.log(fpkm+options.pseudocount,2) for fpkm in gene_fpkm1[cond_key]])
        fpkm2_r = ro.FloatVector([math.log(fpkm+options.pseudocount,2) for fpkm in gene_fpkm2[cond_key]])
        qvals_r = ro.FloatVector([math.log(qval+1e-15,10) for qval in gene_qvals[cond_key]])

        df = ro.DataFrame({'fpkm1':fpkm1_r, 'fpkm2':fpkm2_r, 'qval':qvals_r})

        # plot
        gp = ggplot2.ggplot(df) + \
            ggplot2.aes_string(x='fpkm1', y='fpkm2', colour='qval') + \
            ggplot2.geom_point(size=1.5, alpha=.3) + \
            ggplot2.scale_x_continuous('%s log2 FPKM' % cond1) + \
            ggplot2.scale_y_continuous('%s log2 FPKM' % cond2) + \
            ggplot2.geom_abline(intercept=0, slope=1, linetype=2) + \
            ggplot2.opts(title=plot_title) + \
            ggplot2.theme_bw()

        #ggplot2.scale_colour_gradient(low='red') + \
        
        # save to file
        grdevices.pdf(file='%s_scatter/%s_%s.pdf' % (options.out_dir_pre,cond1,cond2))
        gp.plot()
        grdevices.dev_off()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
