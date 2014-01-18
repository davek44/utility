#!/usr/bin/env python
from optparse import OptionParser
import math, os, subprocess, pdb, shutil, sys
from scipy.stats import spearmanr
import gff, ggplot, stats

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
    parser.add_option('-o', dest='out_dir', default='scatters', help='Prefix for output directories [Default: %default]')
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

        if status == 'OK' and not math.isnan(test_stat) and (len(gene_set) == 0 or gene_id in gene_set):
            gene_qvals.setdefault((cond1,cond2),[]).append(qval)
            gene_fpkm1.setdefault((cond1,cond2),[]).append(fpkm1)
            gene_fpkm2.setdefault((cond1,cond2),[]).append(fpkm2)

        line = diff_in.readline()
    diff_in.close()

    # clean plot directory
    if os.path.isdir('%s_scatter' % options.out_dir_pre):
        shutil.rmtree('%s_scatter' % options.out_dir_pre)
    os.mkdir('%s_scatter' % options.out_dir_pre)

    for cond_key in gene_qvals:
        cond1, cond2 = cond_key

        df_dict = {}
        df_dict['fpkm1'] = [math.log(fpkm+options.pseudocount,2) for fpkm in gene_fpkm1[cond_key]]
        df_dict['fpkm2'] = [math.log(fpkm+options.pseudocount,2) for fpkm in gene_fpkm2[cond_key]]
        df_dict['qval'] = [math.log(qval+1e-15,10) for qval in gene_qvals[cond_key]]

        rho, pval = spearmanr(df_dict['fpkm1'], df_dict['fpkm2'])
        print '%-15s  %-15s  %.4f' % (cond1, cond2, rho)

        output_pdf = '%s/%s_%s.pdf' % (options.out_dir,cond1,cond2)

        ggplot.plot('%s/cuff_scatter.r' % os.environ['RDIR'], df_dict, [output_pdf,cond1,cond2])


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
