#!/usr/bin/env python
from optparse import OptionParser
from numpy import array, empty
from scipy.stats import spearmanr, norm
import random, sys

import rpy2
from rpy2.robjects.numpy2ri import numpy2ri
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

ro.conversion.py2ri = numpy2ri

################################################################################
# array_cors.py
#
# Compute correlations for the genes in a given file, which should use RefSeq
# and lnc_catalog names.
#
# The array data should be in the directory options.cel_dir, which is the 
# current directory by default.
#
# The updated CDF file should be options.cdf, which is 'genes_lnc.cdf' by
# default.
#
# The output will be a file for each gene containing its correlations and then
# another file containing correlations for randomly sampled genes (not from the
# given set).
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <genes file>'
    parser = OptionParser(usage)
    parser.add_option('--cdf', dest='cdf', default='genes_lnc.cdf', help='Chip description file including lncRNAs [Default: %default]')
    parser.add_option('--cel', dest='cel_dir', default='.', help='CEL file directory [Default: %default]')
    parser.add_option('-n', dest='null_samples', type='int', default=200, help='Number of null samples to compute [Default: %default]')
    (options,args) = parser.parse_args()
    
    if len(args) != 1:
        parser.error('Must provide file with names of genes of interest')
    else:
        genes_file = args[0]

    # get genes
    genes_interest = set([line.rstrip() for line in open(genes_file)])

    # determine gene indexes in cdf
    probesets = [line[5:].rstrip() for line in open(options.cdf) if line.startswith('Name=') and line[5:9] != 'NONE']
    probesets = probesets[1:] # first Name= is for the chip
    probesets.sort()
    genes_interest_i = [i for i in range(len(probesets)) if probesets[i][:-3] in genes_interest]

    # choose null gene indexes
    genes_uninterest_i = list(set(range(len(probesets))) - set(genes_interest_i))
    genes_null_i = sorted(random.sample(genes_uninterest_i, options.null_samples))

    # setup R libraries`
    affy = importr('affy')
    makecdfenv = importr('makecdfenv')

    # get lncRNA cdf
    #cdfenv = makecdfenv.make_cdf_env('genes_lnc.cdf')
    cdfenv = ro.r('cdfenv = make.cdf.env("genes_lnc.cdf")')

    # load data
    #array = affy.ReadAffy(cdfname='cdfenv') # or celfile_path="..."
    array = ro.r('array = ReadAffy(celfile.path="%s", cdfname="cdfenv")' % options.cel_dir)

    # run RMA
    #array_rma = affy.rma(array)
    array_rma = ro.r('array_rma = rma(array)')

    # compute correlations
    indexes_str = ','.join([str(ind) for ind in (genes_interest_i+genes_null_i)])
    gene_cors = ro.r('gene_cors = cor(t(exprs(array_rma)), t(exprs(array_rma)[c(%s),]), method="spearman")' % indexes_str)

    # print genes of interest correlations
    for j in range(len(genes_interest_i)):
        i = genes_interest_i[j]
        gint_out = open('%s_cors.txt' % probesets[i][:-3], 'w')
        jcors = gene_cors.rx(True,j+1)
        for ci in range(len(jcors)):
            print >> gint_out, probesets[ci][:-3], jcors[ci]
        gint_out.close()

    # print null genes correlations
    null_out = open('null_cors.txt', 'w')
    for j in range(len(genes_null_i)):
        print >> null_out, '>%d' % j
        jcors = gene_cors.rx(True,j+1+len(genes_interest_i))
        for ci in range(len(jcors)):
            print >> null_out, probesets[ci][:-3], jcors[ci]
    null_out.close()



################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
