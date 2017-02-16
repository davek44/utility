#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
from numpy import array, empty
from scipy.stats import norm
import math, sys
import stats

'''
import rpy2
from rpy2.robjects.numpy2ri import numpy2ri
import rpy2.robjects as ro

ro.conversion.py2ri = numpy2ri
'''

################################################################################
# cufflinks.py
#
# Code to support analysis of cufflinks output.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()



################################################################################
# hash_fpkm
#
# Quick and dirty version to get at a single FPKM value.
################################################################################
def hash_fpkm(fpkm_file, experiment='', fail=float('nan')):
    gene_fpkm = {}

    # get headers
    fpkm_in = open(fpkm_file)
    headers = fpkm_in.readline().split()

    # find experiment column
    exp_col = 0
    while headers[exp_col] != 'FPKM' and headers[exp_col] != '%s_FPKM' % experiment:
        exp_col += 1

    if headers[exp_col] != 'FPKM' and headers[exp_col] != '%s_FPKM' % experiment:
        print >> sys.stderr, '%s unfound' % experiment
        exit(1)

    for line in fpkm_in:
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        gene_id = a[0]

        if a[exp_col+3] in ['FAIL','HIDATA']:
            fpkm = fail
        else:
            fpkm = float(a[exp_col])

        gene_fpkm[gene_id] = fpkm

    fpkm_in.close()

    return gene_fpkm


################################################################################
# hash_fpkms
#
# Quick and dirty version to get the arithmetic mean of a few FPKM values.
################################################################################
def hash_fpkms(fpkm_file, experiments, fail=float('nan')):
    gene_fpkm = {}

    # get headers
    fpkm_in = open(fpkm_file)
    headers = fpkm_in.readline().split()

    # find experiment columns
    exp_cols = []
    for i in range(len(headers)):
        if headers[i][-5:] == '_FPKM':
            if headers[i][:-5] in experiments:
                exp_cols.append(i)

    if len(exp_cols) != len(experiments):
        print >> sys.stderr, '%s unfound' % (','.join(experiments))
        exit(1)

    for line in fpkm_in:
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        gene_id = a[0]

        nonfails = 0
        for exp_col in exp_cols:
            if a[exp_col+3] in ['FAIL','HIDATA']:
                fpkm = fail
            else:
                fpkm = float(a[exp_col])
                nonfails += 1

            gene_fpkm[gene_id] = gene_fpkm.get(gene_id,0) + fpkm

        if nonfails > 0:
            gene_fpkm[gene_id] /= float(nonfails)
        else:
            gene_fpkm[gene_id] = fail

    fpkm_in.close()

    return gene_fpkm


################################################################################
# fpkm_tracking
################################################################################
class fpkm_tracking:
    ############################################################################
    # Constructor
    #
    # Load the expression matrix 
    ############################################################################
    def __init__(self, fpkm_file):
        # obtain basic information
        fpkm_in = open(fpkm_file)
        headers = fpkm_in.readline().split()
        self.genes = []
        line = fpkm_in.readline()
        while line:
            a = line.split()
            self.genes.append(a[0])
            line = fpkm_in.readline()
        fpkm_in.close()

        self.gene_map = dict([(self.genes[i],i) for i in range(len(self.genes))])
        self.experiments = [h[:-5] for h in headers if h == 'FPKM' or h[-5:] == '_FPKM']
        if len(self.experiments) == 1 and self.experiments[0] == '':
            self.experiments[0] = 'unknown'

        # obtain expression
        self.expr = empty([len(self.genes), len(self.experiments)])
        g = 0

        fpkm_in = open(fpkm_file)
        line = fpkm_in.readline()
        line = fpkm_in.readline()
        while line:
            a = line.split('\t')
            a[-1] = a[-1].rstrip()
            e = 0
            for i in range(len(headers)):
                if headers[i] == 'FPKM' or headers[i][-5:] == '_FPKM':
                    if a[i+3] in ['FAIL','HIDATA']:
                        self.expr[g,e] = float('nan')
                    else:
                        self.expr[g,e] = float(a[i])
                    e += 1
            g += 1
            line = fpkm_in.readline()
        fpkm_in.close()

        print >> sys.stderr, 'Loaded expression of %d genes in %d experiments' % (g,e)


    ############################################################################
    # gene_entropy
    #
    # Return the entropy of the expression vector for the given gene.
    # Note:
    #  -Log creates negative values so it's not a distribution for sure.
    ############################################################################
    def gene_entropy(self, gene, log=False):
        gene_i = self.name_or_index(gene)

        gexpr = self.expr[gene_i,:]
        if log:
            gexpr = [math.log(e+1) for e in gexpr]
        gexpr = stats.normalize(gexpr)

        return stats.entropy(gexpr)


    ############################################################################
    # gene_expr
    #
    # Return an expression vector for the given gene.
    ############################################################################
    def gene_expr(self, gene, not_found=float('nan'), fail=float('nan')):
        gene_i = self.name_or_index(gene)
        if gene_i:
            expr_vec = [e if not math.isnan(e) else fail for e in self.expr[gene_i,:]]
            return expr_vec
        else:
            return [not_found]*len(self.experiments)


    ############################################################################
    # gene_expr_exp
    #
    # Return FPKM for a given gene in a given experiment.
    ############################################################################
    def gene_expr_exp(self, gene, exp, not_found=float('nan'), fail=float('nan')):
        gene_i = self.name_or_index(gene)
        if gene_i == None:
            return not_found
        else:
            for exp_i in range(len(self.experiments)):
                if self.experiments[exp_i] == exp:
                    if math.isnan(self.expr[gene_i,exp_i]):
                        return fail
                    else:
                        return self.expr[gene_i,exp_i]
            return not_found


    ############################################################################
    # gene_expr_print
    #
    # Print expression data for the given gene.
    ############################################################################
    def gene_expr_print(self, gene):
        gene_i = self.name_or_index(gene)
            
        for j in range(len(self.experiments)):
            print('%-15s %8.3f' % (self.experiments[j], self.expr[gene_i,j]))


    ############################################################################
    # gene_specificity
    #
    # Return tissue specificity for the given gene.
    ############################################################################
    def gene_specificity(self, gene, log=True):
        gene_i = self.name_or_index(gene)

        if gene_i == None:
            spec = 0
        else:
            gexpr = self.expr[gene_i,:]
            if log:
                gexpr = [math.log(e+1) for e in gexpr]

            if sum(gexpr) == 0:
                spec = 0
            else:
                gexpr = stats.normalize(gexpr)

                min_jsd = 1.0
                for j in range(len(self.experiments)):
                    q_j = [0]*len(self.experiments)
                    q_j[j] = 1.0

                    min_jsd = min(min_jsd, math.sqrt(stats.jsd(gexpr, q_j)))

                spec = 1.0 - min_jsd

        return spec


    ############################################################################
    # genes_jsd
    #
    # Jensen-Shannon divergence between two genes
    ############################################################################
    def genes_jsd(self, gene1, gene2, log=True):
        gene1_i = self.name_or_index(gene1)
        gene2_i = self.name_or_index(gene2)

        gexpr1 = self.expr[gene1_i,:]
        if log:
            gexpr1 = [math.log(e+1) for e in gexpr1]
        gexpr1 = stats.normalize(gexpr1)

        gexpr2 = self.expr[gene2_i,:]
        if log:
            gexpr2 = [math.log(e+1) for e in gexpr2]
        gexpr2 = stats.normalize(gexpr2)

        return stats.jsd(gexpr1, gexpr2)
            

    ############################################################################
    # name_or_index
    #
    # Given a name or index, return an index
    ############################################################################
    def name_or_index(self, gene):
        if type(gene) == str:
            if gene in self.gene_map:
                return self.gene_map[gene]
            else:
                #print >> sys.stderr, 'Missing gene - %s' % gene
                return None
        elif type(gene) == int:
            return gene
        else:
            #print >> sys.stderr, 'Bad gene input'
            return None

    
    ############################################################################
    # spearman
    #
    # Compute Spearman correlations for either all pairs of genes or one given
    # gene to all others.
    ############################################################################
    '''
    def spearman(self, gene=None):
        if type(gene) == str:
            cors = ro.r.cor(self.expr[self.gene_map[gene],:], self.expr.transpose(), method='spearman')
            return array(cors)[0]

        elif type(gene) == int:
            cors = ro.r.cor(self.expr[gene,:], self.expr.transpose(), method='spearman')
            return array(cors)[0]

        else:
            cors = ro.r.cor(self.expr.transpose(), method='spearman')
            return array(cors)
    '''

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
