#!/usr/bin/env python
from optparse import OptionParser
import copy

import numpy as np

'''
quantile_normalization.py
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    x = np.random.randn(10,5)
    for ti in range(x.shape[1]):
        x[:,ti] = (ti+1)*x[:,ti]

    print(x, end='\n\n')

    xn = quantile_normalize(x)
    print(xn, end='\n\n')

    print(x.mean(axis=0))
    print(xn.mean(axis=0))


def quantile_normalize_expr(gene_expr, quantile_stat='median'):
    ''' Quantile normalize across targets. The version below
        just labels the variables more generally, but should
        return the same answer. '''

    # make a copy
    gene_expr_qn = copy.copy(gene_expr)

    # sort values within each column
    for ti in range(gene_expr.shape[1]):
        gene_expr_qn[:,ti].sort()

    # compute the mean/median in each row
    if quantile_stat == 'median':
        sorted_index_stats = np.median(gene_expr_qn, axis=1)
    elif quantile_stat == 'mean':
        sorted_index_stats = np.mean(gene_expr_qn, axis=1)
    else:
        print('Unrecognized quantile statistic %s' % quantile_stat, file=sys.stderr)
        exit()

    # set new values
    for ti in range(gene_expr.shape[1]):
        sorted_indexes = np.argsort(gene_expr[:,ti])
        for gi in range(gene_expr.shape[0]):
            gene_expr_qn[sorted_indexes[gi],ti] = sorted_index_stats[gi]

    return gene_expr_qn

def quantile_normalize(X, quantile_stat='median'):
    ''' Quantile normalize features across samples. '''

    # make a copy
    Xq = copy.copy(X)

    # sort values within each column
    for fi in range(X.shape[1]):
        Xq[:,fi].sort()

    # compute the mean/median in each row
    if quantile_stat == 'median':
        sorted_index_stats = np.median(Xq, axis=1)
    elif quantile_stat == 'mean':
        sorted_index_stats = np.mean(Xq, axis=1)
    else:
        print('Unrecognized quantile statistic %s' % quantile_stat, file=sys.stderr)
        exit()

    # set new values
    for fi in range(X.shape[1]):
        sorted_indexes = np.argsort(X[:,fi])
        for si in range(X.shape[0]):
            Xq[sorted_indexes[si],fi] = sorted_index_stats[si]

    return Xq


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
