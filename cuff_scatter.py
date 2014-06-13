#!/usr/bin/env python
from optparse import OptionParser
import math, os, subprocess, pdb, shutil, sys
from scipy.stats import spearmanr, pearsonr
import gff, ggplot, stats
import cufflinks

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
    parser.add_option('-f', dest='fpkm_tracking', default=False, action='store_true', help='Providing an fpkm_tracking file rather than a diff file [Default: %default]')
    parser.add_option('-g', dest='genes_gtf', help='Print only genes in the given GTF file')
    parser.add_option('-o', dest='out_dir', default='scatters', help='Prefix for output directories [Default: %default]')
    parser.add_option('-p', dest='pseudocount', type='float', default=0.125, help='FPKM pseudocount for taking logs [Default: %default]')
    parser.add_option('-r', dest='read_group_tracking', default=False, action='store_true', help='Providing a reads_group_tracking file rather than a diff file [Default: %default]')
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

    conditions_gene_qval = {}
    if options.fpkm_tracking:
        condition_gene_fpkm = read_fpkm_tracking(diff_file, gene_set)
    elif options.read_group_tracking:
        condition_gene_fpkm = read_read_group_tracking(diff_file, gene_set)
    else:
        condition_gene_fpkm, conditions_gene_qval = read_diff(diff_file, gene_set)

    # fill in dummy q-value conditions
    if len(conditions_gene_qval) == 0:
        conditions = condition_gene_fpkm.keys()
        for i in range(len(conditions)):
            for j in range(i+1,len(conditions)):
                if options.read_group_tracking:
                    cond1_pre = conditions[i][:conditions[i].rfind('rep')]
                    cond2_pre = conditions[j][:conditions[j].rfind('rep')]
                    if cond1_pre == cond2_pre:
                        conditions_gene_qval[(conditions[i],conditions[j])] = {}
                else:
                    conditions_gene_qval[(conditions[i],conditions[j])] = {}

    # clean plot directory
    if os.path.isdir(options.out_dir):
        shutil.rmtree(options.out_dir)
    os.mkdir(options.out_dir)

    for cond_key in conditions_gene_qval:
        cond1, cond2 = cond_key

        if len(conditions_gene_qval[cond_key]) > 0:
            plot_genes = conditions_gene_qval[cond_key].keys()
        else:
            plot_genes = set(condition_gene_fpkm[cond1].keys()) & set(condition_gene_fpkm[cond2].keys())

        df_dict = {'fpkm1': [math.log(condition_gene_fpkm[cond1][gene_id]+options.pseudocount,2) for gene_id in plot_genes],
                   'fpkm2': [math.log(condition_gene_fpkm[cond2][gene_id]+options.pseudocount,2) for gene_id in plot_genes],
                   'qval': [math.log(conditions_gene_qval[cond_key].get(gene_id,1-1e-15)+1e-15,10) for gene_id in plot_genes]}        

        #rho, pval = spearmanr(df_dict['fpkm1'], df_dict['fpkm2'])
        rho, pval = pearsonr(df_dict['fpkm1'], df_dict['fpkm2'])
        print '%-15s  %-15s  %.4f' % (cond1, cond2, rho)

        output_pdf = '%s/%s_%s.pdf' % (options.out_dir,cond1,cond2)

        ggplot.plot('%s/cuff_scatter.r' % os.environ['RDIR'], df_dict, [output_pdf,cond1,cond2,options.pseudocount])


################################################################################
# read_diff
################################################################################
def read_diff(diff_file, gene_set):
    # initialize diff data structures
    condition_gene_fpkm = {}
    conditions_gene_qval = {}

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

        if cond2 in ['input','control']:
            cond1, cond2 = cond2, cond1
            fpkm1, fpkm2 = fpkm2, fpkm1
            fold_change *= -1
            test_stat *= -1
            
        if status == 'OK' and not math.isnan(test_stat) and (len(gene_set) == 0 or gene_id in gene_set):
            conditions_gene_qval.setdefault((cond1,cond2),{})[gene_id] = qval
            condition_gene_fpkm.setdefault(cond1,{})[gene_id] = fpkm1
            condition_gene_fpkm.setdefault(cond2,{})[gene_id] = fpkm2

        line = diff_in.readline()
    diff_in.close()

    return condition_gene_fpkm, conditions_gene_qval


################################################################################
# read_fpkm_tracking
################################################################################
def read_fpkm_tracking(ft_file, gene_set):
    condition_gene_fpkm = {}

    # get fpkm's
    cuff = cufflinks.fpkm_tracking(ft_file)

    # for each gene
    for gene_i in range(len(cuff.genes)):
        gene_id = cuff.genes[gene_i]
        if len(gene_set) == 0 or gene_id in gene_set:
            gene_expr = cuff.gene_expr(gene_i)

            # for each condition
            for i in range(len(cuff.experiments)):
                # save, as long as not nan
                if not math.isnan(gene_expr[i]):
                    cond = cuff.experiments[i]
                    condition_gene_fpkm.setdefault(cond,{})[gene_id] = gene_expr[i]

    return condition_gene_fpkm


################################################################################
# read_read_group_tracking
################################################################################
def read_read_group_tracking(rgt_file, gene_set):
    condition_gene_fpkm = {}

    rgt_in = open(rgt_file)
    line = rgt_in.readline()
    for line in rgt_in:
        a = line.split('\t')
        gene_id = a[0]
        condition = a[1]
        rep = int(a[2])
        fpkm = float(a[6])

        if len(gene_set) == 0 or gene_id in gene_set:
            cond_key = '%s-rep%d' % (condition,rep)
            condition_gene_fpkm.setdefault(cond_key,{})[gene_id] = fpkm            

    return condition_gene_fpkm

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
