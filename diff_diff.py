#!/usr/bin/env python
from optparse import OptionParser
import os, math, sys
from scipy.stats import spearmanr
import ggplot, ripseq

################################################################################
# diff_diff.py
#
# Compare two cuffdiff runs.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <diff1_file> <diff2_file>'
    parser = OptionParser(usage)
    parser.add_option('-o', dest='out_dir', default='.')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide two diff files')
    else:
        diff1_file = args[0]
        diff2_file = args[1]

    diff1_stats, diff1_bound = hash_diff(diff1_file)
    diff2_stats, diff2_bound = hash_diff(diff2_file)

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    for diff_key in diff1_stats:
        sample1, sample2 = diff_key

        gene_stats1 = diff1_stats[diff_key]
        gene_bound1 = diff1_bound[diff_key]
        gene_stats2 = diff2_stats[diff_key]
        gene_bound2 = diff2_bound[diff_key]

        report_out = open('%s/%s-%s_report.txt' % (options.out_dir,sample1,sample2), 'w')

        # compare numbers of genes quantified
        common_genes = set(gene_stats1.keys()) & set(gene_stats2.keys())
        print >> report_out, 'Genes quantified'
        print >> report_out, '%s\t%d' % (diff1_file,len(gene_stats1))
        print >> report_out, '%s\t%d' % (diff2_file,len(gene_stats2))
        print >> report_out, 'Common\t%d' % len(common_genes)
        print >> report_out, ''

        up1 = set([gene_id for gene_id in gene_bound1 if gene_bound1[gene_id]])
        up2 = set([gene_id for gene_id in gene_bound2 if gene_bound2[gene_id]])

        print >> report_out, 'Genes upregulated'    
        print >> report_out, '%s\t%d' % (diff1_file,len(up1))
        print >> report_out, '%s\t%d' % (diff2_file,len(up2))
        print >> report_out, 'Common\t%d' % len(up1 & up2)
        print >> report_out, ''

        down1 = set([gene_id for gene_id in gene_bound1 if not gene_bound1[gene_id]])
        down2 = set([gene_id for gene_id in gene_bound2 if not gene_bound2[gene_id]])

        print >> report_out, 'Genes downregulated'    
        print >> report_out, '%s\t%d' % (diff1_file,len(down1))
        print >> report_out, '%s\t%d' % (diff2_file,len(down2))
        print >> report_out, 'Common\t%d' % len(down1 & down2)
        print >> report_out, ''

        # scatter plot test stat
        df = {'diff1':[], 'diff2':[]}
        for gene_id in common_genes:
            df['diff1'].append(gene_stats1[gene_id])
            df['diff2'].append(gene_stats2[gene_id])

        r_script = '%s/diff_diff_scatter.r' % os.environ['RDIR']
        out_pdf = '%s/%s-%s_scatter.pdf' % (options.out_dir, sample1, sample2)
        ggplot.plot(r_script, df, [out_pdf])

        # compute correlation
        cor, p = spearmanr(df['diff1'], df['diff2'])

        print >> report_out, 'Spearman correlation: %f' % cor
        print >> report_out, ''

        report_out.close()

        # plot test_stat versus test_stat difference
        df = {'minus':[], 'avg':[]}
        for gene_id in common_genes:
            df['minus'].append(gene_stats1[gene_id] - gene_stats2[gene_id])
            df['avg'].append(0.5*gene_stats1[gene_id] + 0.5*gene_stats2[gene_id])

        r_script = '%s/diff_diff_ma.r' % os.environ['RDIR']
        out_pdf = '%s/%s-%s_ma.pdf' % (options.out_dir, sample1, sample2)
        ggplot.plot(r_script, df, [out_pdf])


################################################################################
# hash_diff
# 
# Input:
#  diff_file:   gene_exp.diff file
#  min_fpkm:    Minimum FPKM to consider a gene.
#
# Output:
#  gene_tstat:  Dict mapping sample pairs to dicts mapping gene_id to test_stat
#  gene_bound:  Dict mapping sample pairs to dicts mapping gene_id to True/False
################################################################################
def hash_diff(diff_file, use_fold=False, min_fpkm=None):
    diff_stat = {}
    diff_bound = {}

    # read rip diff
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
        fold_change = float(a[9])
        tstat = float(a[10])
        qval = float(a[11])
        sig = a[-1].rstrip()

        if sample1 > sample2:
            sample1, sample2 = sample2, sample1
            fpkm1, fpkm2 = fpkm2, fpkm1
            fold_change *= -1
            tstat *= -1

        diff_key = (sample1, sample2)

        if status == 'OK' and not math.isnan(tstat):
            if min_fpkm == None or fpkm1 > min_fpkm or fpkm2 > min_fpkm:
                if use_fold:
                    diff_stat.setdefault(diff_key,{})[gene_id] = fold_change
                else:
                    diff_stat.setdefault(diff_key,{})[gene_id] = tstat

                if sig == 'yes':
                    if tstat > 0:
                        diff_bound.setdefault(diff_key,{})[gene_id] = True
                    else:
                        diff_bound.setdefault(diff_key,{})[gene_id] = False

    diff_in.close()

    # add sample pairs w/ no sig genes
    for diff_key in diff_stat:
        if not diff_key in diff_bound:
            diff_bound[diff_key] = {}

    return diff_stat, diff_bound


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
