#!/usr/bin/env python
from optparse import OptionParser
import copy, math, os, pdb, random, sys

import numpy as np
from sklearn import preprocessing
from sklearn.decomposition import PCA, FastICA
from sklearn.manifold import MDS, Isomap

import gff, ggplot

################################################################################
# cuff_2d.py
#
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <read_group_tracking>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='gtf', help='GTF file of genes to display')
    parser.add_option('-f', dest='min_fpkm', default=0, type='float', help='Minimum FPKM to consider [Default: %default]')
    parser.add_option('-m', dest='method', default='PCA', help='Dimension reduction method [Default: %default]')
    parser.add_option('-p', dest='pseudocount', default=.125, help='FPKM pseudocount (for logs) [Default: %default]')
    parser.add_option('-o', dest='out_pdf', default='cuff_2d.pdf', help='Output PDF [Default: %default]')
    parser.add_option('-s', dest='square', default=False, action='store_true', help='Square plot [Default: %default]')
    parser.add_option('-u', dest='uppercase', default=False, action='store_true', help='Uppercase sample labels [Default: %default]')
    parser.add_option('-w', dest='whiten', default=False, action='store_true', help='Whiten expression data [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide fpkm_tracking')
    else:
        read_group_tracking = args[0]

    # load expression data
    gene_fpkm = {}
    rgt_in = open(read_group_tracking)
    rgt_in.readline()
    for line in rgt_in:
        a = line.split()
        gene_fpkm.setdefault(a[0],{})[(a[1],a[2])] = float(a[6])
    rgt_in.close()

    # determine genes
    compute_genes = gene_fpkm.keys()
    if options.gtf:
        compute_genes = set()
        for line in open(options.gtf):
            a = line.split('\t')
            compute_genes.add(gff.gtf_kv(a[8])['gene_id'])
    compute_genes = list(compute_genes)

    # filter for fpkm
    if options.min_fpkm > 0:
        prefilter_genes = copy.copy(compute_genes)
        compute_genes = []
        for gene_id in prefilter_genes:
            ge = gene_fpkm[gene_id].values()
            if max(ge) > options.min_fpkm:
                compute_genes.append(gene_id)

    # construct gene expression matrix
    samples = gene_fpkm[compute_genes[0]].keys()
    X = np.array([[gene_fpkm[gene_id][sam_rep] for gene_id in compute_genes] for sam_rep in samples])
    X = np.log2(X + options.pseudocount)

    if options.whiten:
        X = preprocessing.scale(X)

    ##################################################
    # dimensionality reduction
    ##################################################
    if options.method.lower() == 'mds':
        model = MDS(n_components=2)
    elif options.method.lower() in ['iso','isomap']:
        model = Isomap(n_components=2)
    elif options.method.lower() == 'ica':
        model = FastICA(n_components=2, max_iter=500)
    else:
        model = PCA(n_components=2)
    
    X_dr = model.fit_transform(X)

    ##################################################    
    # plot
    ##################################################
    df = {}
    df['D1'] = X_dr[:,0]
    df['D2'] = X_dr[:,1]
    df['Label'] = ['%s_rep%s' % sam_rep for sam_rep in samples]
     
    if options.uppercase:
        df['Label'] = [label.upper() for label in df['Label']]
        df['Sample'] = [sam.upper() for (sam,rep) in samples]
    else:
        df['Sample'] = [sam for (sam,rep) in samples]

    ggplot.plot('%s/cuff_2d.r' % os.environ['RDIR'], df, [options.out_pdf, options.square])


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
