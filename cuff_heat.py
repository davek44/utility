#!/usr/bin/env python
from optparse import OptionParser
import numpy as np
import math, os, pdb, random, sys
import cufflinks, gff, ggplot

################################################################################
# cuff_heat.py
#
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fpkm_tracking>'
    parser = OptionParser(usage)
    parser.add_option('-a', dest='max_fpkm', type='float', help='Maxium log2 FPKM to plot [Default: %d]')
    parser.add_option('-d', dest='diff_file', help='Limit to significantly differentially expressed genes')
    parser.add_option('-g', dest='gtf', help='GTF file of genes to display')
    parser.add_option('-m', dest='min_fpkm', default=0, type='float', help='Minimum FPKM [Default: %default]')
    parser.add_option('-p', dest='pseudocount', default=.125, type='float', help='Pseudocount for log FPKM [Default: %default]')
    parser.add_option('-o', dest='out_pdf', default='cuff_heat.pdf', help='Output PDF [Default: %default]')
    parser.add_option('-s', dest='sample', default=1000, type='int', help='Sample genes rather than use all [Default: %default]')
    parser.add_option('-u', dest='uppercase', default=False, action='store_true', help='Uppercase sample labels [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide fpkm_tracking')
    else:
        fpkm_tracking = args[0]

    # load expression data
    cuff = cufflinks.fpkm_tracking(fpkm_file=fpkm_tracking)

    # determine genes
    all_genes = set(cuff.genes)
    if options.gtf:
        all_genes = set()
        for line in open(options.gtf):
            a = line.split('\t')
            all_genes.add(gff.gtf_kv(a[8])['gene_id'])


    if options.diff_file:
        # limit to differentially expressed genes
        diff_genes = find_diff(options.diff_file)
        all_genes &= diff_genes

    else:
        # at least limit to clean genes
        clean_genes = set()
        for gene_id in all_genes:
            ge = cuff.gene_expr(gene_id)
            clean = True 
            for i in range(len(ge)):
                if math.isnan(ge[i]):
                    clean = False
                    break
            if clean:
                clean_genes.add(gene_id)

        all_genes &= clean_genes


    if options.min_fpkm > 0:
        expressed_genes = set()
        for gene_id in all_genes:
            ge = cuff.gene_expr(gene_id, not_found=0, fail=0)
            if max(ge) > options.min_fpkm:
                expressed_genes.add(gene_id)
                
        all_genes &= expressed_genes

    # sample genes to display
    if len(all_genes) <= options.sample:
        display_genes = all_genes
    else:
        display_genes = random.sample(all_genes, options.sample)

    # build data frame
    df = {'Gene':[], 'FPKM':[], 'Sample':[]}

    for gene_id in display_genes:
        ge = cuff.gene_expr(gene_id, not_found=0, fail=0)

        for i in range(len(cuff.experiments)):
            df['Gene'].append(gene_id)

            df['Sample'].append(cuff.experiments[i])
            if options.uppercase:
                df['Sample'][-1] = df['Sample'][-1].upper()

            logfpkm = np.log2(ge[i]+options.pseudocount)
            if options.max_fpkm:
                logfpkm = min(options.max_fpkm, logfpkm)
            df['FPKM'].append(logfpkm)

    # plot
    out_df = '%s.df' % options.out_pdf[:-4]
    ggplot.plot('%s/cuff_heat.r' % os.environ['RDIR'], df, [options.out_pdf], df_file=out_df)


################################################################################
# find_diff
#
# Return a set of only the differentially expressed genes.
################################################################################
def find_diff(diff_file):
    diff_genes = set()

    diff_in = open(diff_file)
    diff_in.readline()

    for line in diff_in:
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        gene_id = a[0]
        sig = a[-1]

        if sig == 'yes':
            diff_genes.add(gene_id)

    diff_in.close()

    return diff_genes


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
