#!/usr/bin/env python
from optparse import OptionParser
import math, os, pdb, shutil, subprocess

import numpy as np
import pandas as pd

import cuffdiff, fdr, gff, ggplot, te, stats

################################################################################
# te_diff_regress.py
#
# Run a regression to determine TE affect on differential expression.
################################################################################

count_tes = ['LINE/L1', 'LINE/L2', 'SINE/Alu', 'SINE/MIR', 'LTR/ERV1', 'LTR/ERVL', 'LTR/ERVL-MaLR', 'DNA/hAT-Charlie', 'DNA/TcMar-Tigger']
#regression_tes = ['LINE/L1', 'LINE/L2', 'LINE/CR1', 'LINE/RTE-X', 'SINE/Alu', 'SINE/MIR', 'LTR/ERV1', 'LTR/ERVL', 'LTR/ERVL-MaLR', 'LTR/Gypsy', 'LTR/ERVK', 'DNA/hAT-Charlie', 'DNA/TcMar-Tigger', 'DNA/hAT-Tip100', 'DNA/hAT-Blackjack', 'DNA/TcMar-Mariner']


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf> <diff>'
    parser = OptionParser(usage)
    parser.add_option('-o', dest='out_dir', default='te_diff_regress', help='Output directory to print regression summaries [Default: %default]')
    parser.add_option('-c', dest='scale', default=1, type='float', help='Plot scale [Default: %default]')


    parser.add_option('-t', dest='te_gff', default='%s/hg19.fa.out.tp.gff'%os.environ['MASK'])
    parser.add_option('-r', dest='orientation', default=False, action='store_true', help='Split TEs by orientation [Default: %default]')

    parser.add_option('-m', dest='max_stat', default=None, type='float', help='Maximum stat for plotting [Default: %default]')

    parser.add_option('-s', dest='spread_factor', default=None, type='float', help='Allow multiplicative factor between the shortest and longest transcripts, used to filter [Default: %default]')
    parser.add_option('-l', dest='spread_lower', default=None, type='float', help='Allow multiplicative factor between median and shortest transcripts [Defafult: %default]')
    parser.add_option('-u', dest='spread_upper', default=None, type='float', help='Allow multiplicative factor between median and longest transcripts [Defafult: %default]')

    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide .gtf and .diff files')
    else:
        ref_gtf = args[0]
        diff_file = args[1]

    # make output directory
    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    if options.spread_factor or options.spread_lower or options.spread_upper:
        # filter for similar length

        if options.spread_factor:
            options.spread_lower = math.sqrt(options.spread_factor)
            options.spread_upper = math.sqrt(options.spread_factor)

        spread_gtf = '%s/spread_factor.gtf' % options.out_dir
        gff.length_filter(ref_gtf, spread_gtf, options.spread_lower, options.spread_lower, verbose=True)

        ref_gtf = spread_gtf

    # hash genes -> TEs -> occurence num
    gene_te_num = te.hash_genes_repeats_num(ref_gtf, options.te_gff, gene_key='transcript_id', add_star=True, stranded=options.orientation)

    # hash diffs stats
    gene_diffs = cuffdiff.hash_diff(diff_file, stat='fold', max_stat=options.max_stat, sample_first='input')

    table_lines = []
    pvals = []

    for spair in gene_diffs:
        sample1, sample2 = spair

        gene_list = list(set(gene_te_num.keys()) & set(gene_diffs[spair].keys()))
                
        for fam in count_tes:
            if options.orientation:
                orients = ['+','-']
            else:
                orients = ['+']

            for orient in orients:
                # hash diff values by TE count
                count_diff = []
                for gene_id in gene_diffs[spair]:
                    if options.orientation:
                        count = gene_te_num.get(gene_id,{}).get(('*',fam,orient), 0)
                    else:
                        count = gene_te_num.get(gene_id,{}).get(('*',fam), 0)

                    while count >= len(count_diff):
                        count_diff.append([])
                    count_diff[count].append(gene_diffs[spair][gene_id])

                df = {'TEs':[], 'stat_low':[], 'stat_mid':[], 'stat_hi':[]}
                for c in range(len(count_diff)):
                    if len(count_diff[c]) > 12:
                        stat_low, stat_mid, stat_hi = stats.quantile(count_diff[c], [.25, .5, .75])
                        df['TEs'].append(c)
                        df['stat_low'].append(stat_low)
                        df['stat_mid'].append(stat_mid)
                        df['stat_hi'].append(stat_hi)
                    else:
                        break

                if len(df['TEs']) > 1:
                    fam_plot = fam[fam.find('/')+1:]

                    if options.orientation:                        
                        out_pdf = '%s/%s-%s_%s_%s.pdf' % (options.out_dir, sample1, sample2, fam_plot, orient)
                        out_df = '%s/%s-%s_%s_%s.df' % (options.out_dir, sample1, sample2, fam_plot, orient)
                    else:
                        out_pdf = '%s/%s-%s_%s.pdf' % (options.out_dir, sample1, sample2, fam_plot)
                        out_df = '%s/%s-%s_%s.df' % (options.out_dir, sample1, sample2, fam_plot)

                    ggplot.plot('%s/te_diff_count.r' % os.environ['RDIR'], df, [out_pdf, options.scale], df_file=out_df)

    if options.spread_factor or options.spread_lower or options.spread_upper:
        os.remove(spread_gtf)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
