#!/usr/bin/env python
from optparse import OptionParser
import os, math, sys
from scipy.stats import spearmanr, pearsonr
import cuffdiff, gff, ggplot, ripseq

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
    parser.add_option('-s', dest='stat', default='test_stat')
    parser.add_option('-g', dest='genes_gtf', default=None)
    parser.add_option('-m', dest='min_fpkm', default=0, type='float')    
    parser.add_option('-o', dest='out_dir', default='.')
    parser.add_option('-p', dest='pseudocount', default=0.125, type='float')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide two diff files')
    else:
        diff1_file = args[0]
        diff2_file = args[1]

    gtf_genes = None
    if options.genes_gtf:
        gtf_genes = gff.gtf_gene_set(options.genes_gtf)

    diff1_stats = cuffdiff.hash_stat(diff1_file, stat=options.stat, min_fpkm=options.min_fpkm, pseudocount=options.pseudocount, gene_set=gtf_genes)
    diff1_sig = cuffdiff.hash_sig(diff1_file, gene_set=gtf_genes)

    diff2_stats = cuffdiff.hash_stat(diff2_file, stat=options.stat, min_fpkm=options.min_fpkm, pseudocount=options.pseudocount, gene_set=gtf_genes)
    diff2_sig = cuffdiff.hash_sig(diff2_file, gene_set=gtf_genes)
    
    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    for diff_key in diff1_stats:
        sample1, sample2 = diff_key

        gene_stats1 = diff1_stats[diff_key]
        gene_sig1 = diff1_sig[diff_key]
        gene_stats2 = diff2_stats[diff_key]
        gene_sig2 = diff2_sig[diff_key]

        report_out = open('%s/%s-%s_report.txt' % (options.out_dir,sample1,sample2), 'w')

        # compare numbers of genes quantified
        common_genes = set(gene_stats1.keys()) & set(gene_stats2.keys())
        print >> report_out, 'Genes quantified'
        print >> report_out, '%s\t%d' % (diff1_file,len(gene_stats1))
        print >> report_out, '%s\t%d' % (diff2_file,len(gene_stats2))
        print >> report_out, 'Common\t%d' % len(common_genes)
        print >> report_out, ''

        up1 = set([gene_id for gene_id in gene_sig1 if gene_sig1[gene_id]])
        up2 = set([gene_id for gene_id in gene_sig2 if gene_sig2[gene_id]])

        print >> report_out, 'Genes upregulated'    
        print >> report_out, '%s\t%d' % (diff1_file,len(up1))
        print >> report_out, '%s\t%d' % (diff2_file,len(up2))
        print >> report_out, 'Common\t%d' % len(up1 & up2)
        print >> report_out, ''

        down1 = set([gene_id for gene_id in gene_sig1 if not gene_sig1[gene_id]])
        down2 = set([gene_id for gene_id in gene_sig2 if not gene_sig2[gene_id]])

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
        ggplot.plot(r_script, df, [out_pdf], df_file='%s.df'%out_pdf[:-4])

        # compute correlation
        cor, p = spearmanr(df['diff1'], df['diff2'])
        print >> report_out, 'Spearman correlation: %f (%f)' % (cor,p)
        cor, p = pearsonr(df['diff1'], df['diff2'])
        print >> report_out, 'Pearson correlation: %f (%f)' % (cor,p)

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
# __main__
################################################################################
if __name__ == '__main__':
    main()
