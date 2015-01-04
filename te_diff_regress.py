#!/usr/bin/env python
from optparse import OptionParser
import math, os, pdb, shutil, subprocess

import numpy as np
import pandas as pd

import statsmodels.formula.api as smf

import cuffdiff, fdr, gff, te

################################################################################
# te_diff_regress.py
#
# Run a regression to determine TE affect on differential expression.
################################################################################

regression_tes = ['LINE/L1', 'LINE/L2', 'SINE/Alu', 'SINE/MIR', 'LTR/ERV1', 'LTR/ERVL', 'LTR/ERVL-MaLR', 'DNA/hAT-Charlie', 'DNA/TcMar-Tigger']
#regression_tes = ['LINE/L1', 'LINE/L2', 'LINE/CR1', 'LINE/RTE-X', 'SINE/Alu', 'SINE/MIR', 'LTR/ERV1', 'LTR/ERVL', 'LTR/ERVL-MaLR', 'LTR/Gypsy', 'LTR/ERVK', 'DNA/hAT-Charlie', 'DNA/TcMar-Tigger', 'DNA/hAT-Tip100', 'DNA/hAT-Blackjack', 'DNA/TcMar-Mariner']


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf> <diff>'
    parser = OptionParser(usage)
    parser.add_option('-o', dest='out_dir', default='te_diff_regress', help='Output directory to print regression summaries [Default: %default]')
    parser.add_option('-t', dest='te_gff', default='%s/hg19.fa.out.tp.gff'%os.environ['MASK'])
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

    # clean output directory
    if os.path.isdir(options.out_dir):
        shutil.rmtree(options.out_dir)
    os.mkdir(options.out_dir)

    if options.spread_factor or options.spread_lower or options.spread_upper:
        # filter for similar length

        if options.spread_factor:
            options.spread_lower = math.sqrt(options.spread_factor)
            options.spread_upper = math.sqrt(options.spread_factor)

        spread_gtf = '%s/spread_factor.gtf' % options.out_dir
        gff.length_filter(ref_gtf, spread_gtf, options.spread_lower, options.spread_lower, verbose=True)

        ref_gtf = spread_gtf

    # hash genes -> TEs
    gene_tes = te.hash_genes_repeats(ref_gtf, options.te_gff, gene_key='transcript_id', add_star=True, stranded=True)

    # hash diffs stats
    gene_diffs = cuffdiff.hash_diff(diff_file, stat='fold', max_stat=options.max_stat, sample_first='input')

    table_lines = []
    pvals = []

    for spair in gene_diffs:
        sample1, sample2 = spair

        # construct data frame
        gene_list = list(set(gene_tes.keys()) & set(gene_diffs[spair].keys()))
        df = pd.DataFrame({'diff': [gene_diffs[spair][gene_id] for gene_id in gene_list]})

        covariate_str = ''
        for fam in regression_tes:
            te_key = '%s_fwd' % fam.replace('/','_').replace('-','')
            df[te_key] = [1.0*(('*',fam,'+') in gene_tes.get(gene_id,[])) for gene_id in gene_list]
            if len(covariate_str) == 0:
                covariate_str = te_key
            else:
                covariate_str += ' + %s' % te_key

            te_key = '%s_rev' % fam.replace('/','_').replace('-','')
            df[te_key] = [1.0*(('*',fam,'-') in gene_tes.get(gene_id,[])) for gene_id in gene_list]
            covariate_str += ' + %s' % te_key    

        # regress
        mod = smf.ols(formula='diff ~ %s' % covariate_str, data=df).fit()

        # output model
        mod_out = open('%s/%s-%s.txt' % (options.out_dir, sample1, sample2), 'w')
        print >> mod_out, mod.summary()
        mod_out.close()

        # save table lines
        for fam in regression_tes:
            te_key = '%s_fwd' % fam.replace('/','_').replace('-','')
            cols = (fam, '+', sample1, sample2, sum(df[te_key]), mod.params[te_key], mod.tvalues[te_key], mod.pvalues[te_key]/0.5)
            table_lines.append('%-17s  %1s  %-10s  %-10s  %6d  %8.3f  %8.3f  %10.2e' % cols)
            pvals.append(cols[-1])

            te_key = '%s_rev' % fam.replace('/','_').replace('-','')
            cols = (fam, '-', sample1, sample2, sum(df[te_key]), mod.params[te_key], mod.tvalues[te_key], mod.pvalues[te_key]/0.5)
            table_lines.append('%-17s  %1s  %-10s  %-10s  %6d  %8.3f  %8.3f  %10.2e' % cols)
            pvals.append(cols[-1])

    # perform multiple hypothesis correction
    qvals = fdr.ben_hoch(pvals)

    table_out = open('%s/table.txt' % options.out_dir, 'w')
    for i in range(len(table_lines)):
        print >> table_out, '%s %10.2e' % (table_lines[i],qvals[i])
    table_out.close()


    if options.spread_factor or options.spread_lower or options.spread_upper:
        os.remove(spread_gtf)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
