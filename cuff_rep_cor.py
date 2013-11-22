#!/usr/bin/env python
from optparse import OptionParser
import os, pdb, sys
from numpy import array
from scipy.stats import spearmanr
import gff, ggplot

################################################################################
# cuff_rep_cor.py
#
# Compute correlations between replicates in a cufflinks run.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <.read_group_tracking>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='genes_gtf', help='Print only genes in the given GTF file')
    #parser.add_option('-p', dest='pseudocount', type='float', default=0.125, help='FPKM pseudocount for taking logs [Default: %default]')
    parser.add_option('-o', dest='out_pdf', default='cor_heat.pdf', help='Output heatmap pdf [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error(usage)
    else:
        read_group_tracking = args[0]
    
    # get gene_ids
    gene_set = set()
    if options.genes_gtf:
        for line in open(options.genes_gtf):
            a = line.split('\t')
            gid = gff.gtf_kv(a[8])['gene_id']
            gene_set.add(gid)

    # initialize diff data structures
    cond_rep_gene_fpkm = {}

    # read read group tracking file
    rgt_in = open(read_group_tracking)
    headers = rgt_in.readline()
    line = rgt_in.readline()
    while line:
        a = line.split('\t')

        gene_id = a[0]
        cond = a[1]
        rep = int(a[2])
        fpkm = float(a[6])
        status = a[8].rstrip()

        if status == 'OK' and (len(gene_set) == 0 or gene_id in gene_set):
            if not (cond,rep) in cond_rep_gene_fpkm:
                cond_rep_gene_fpkm[(cond,rep)] = {}
            
            cond_rep_gene_fpkm[(cond,rep)][gene_id] = fpkm

        line = rgt_in.readline()
    rgt_in.close()

    df_dict = {'Sample1':[], 'Sample2':[], 'Correlation':[]}
    cond_reps = cond_rep_gene_fpkm.keys()

    for i in range(len(cond_reps)):
        cond1, rep1 = cond_reps[i]

        for j in range(i+1,len(cond_reps)):
            cond2, rep2 = cond_reps[j]

            genes12 = set(cond_rep_gene_fpkm[(cond1,rep1)].keys()) & set(cond_rep_gene_fpkm[(cond2,rep2)].keys())

            fpkms1 = array([cond_rep_gene_fpkm[(cond1,rep1)][gene_id] for gene_id in genes12])
            fpkms2 = array([cond_rep_gene_fpkm[(cond2,rep2)][gene_id] for gene_id in genes12])

            rho, pval = spearmanr(fpkms1, fpkms2)

            cols = (cond1,rep1,cond2,rep2,rho)
            print '%-15s  %1d  %-15s  %1d  %.4f' % cols

            df_dict['Sample1'].append('%s_%d' % (cond1,rep1))
            df_dict['Sample2'].append('%s_%d' % (cond2,rep2))
            df_dict['Correlation'].append(rho)

    # this is broken
    ggplot.plot('%s/cuff_rep_cor.r' % os.environ['RDIR'], df_dict, [options.out_pdf], debug=True)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
