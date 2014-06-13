#!/usr/bin/env python
from optparse import OptionParser
import math, os, pdb, shutil, subprocess

import numpy as np
import pandas as pd

import statsmodels.formula.api as smf

import fdr, te

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
    parser.add_option('-p', dest='plot_dir', default=None, help='Plot output directory [Default: %default]')
    parser.add_option('-t', dest='te_gff', default='%s/hg19.fa.out.tpf.gff'%os.environ['MASK'])
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide .gtf and .diff files')
    else:
        gtf_file = args[0]
        diff_file = args[1]

    # hash genes -> TEs
    gene_tes = te.hash_genes_repeats(gtf_file, options.te_gff, gene_key='transcript_id', add_star=True, stranded=True)

    # hash diffs stats
    gene_diffs = hash_diff(diff_file)

    # clean plot directory
    if options.plot_dir != None:
        if os.path.isdir(options.plot_dir):
            shutil.rmtree(options.plot_dir)
        os.mkdir(options.plot_dir)

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
        mod_out = open('%s/%s-%s.txt' % (options.plot_dir, sample1, sample2), 'w')
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

    table_out = open('%s/table.txt' % options.plot_dir, 'w')
    for i in range(len(table_lines)):
        print >> table_out, '%s %10.2e' % (table_lines[i],qvals[i])
    table_out.close()


################################################################################
# hash_diff
################################################################################
def hash_diff(diff_file):
    gene_stats = {}

    # read diff file
    diff_in = open(diff_file)
    headers = diff_in.readline()
    line = diff_in.readline()
    while line:
        a = line.split('\t')

        gene_id = a[0]
        sample1 = a[4]
        sample2 = a[5]
        status = a[6]
        fpkm1 = float(a[7])
        fpkm2 = float(a[8])
        fold = float(a[9])
        tstat = float(a[10])
        sig = a[-1].rstrip()

        if sample2 == 'input':
            sample1, sample2 = sample2, sample1
            fpkm1, fpkm2 = fpkm2, fpkm1
            fold *= -1
            tstat *= -1

        # cap fold/tstat
        fold = min(fold, 6)
        fold = max(fold, -6)
        tstat = min(tstat, 6)
        tstat = max(tstat, -6)

        if status == 'OK' and not math.isnan(tstat):
            #gene_stats.setdefault((sample1,sample2),{})[gene_id] = tstat
            gene_stats.setdefault((sample1,sample2),{})[gene_id] = fold

        line = diff_in.readline()
    diff_in.close()

    return gene_stats


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
