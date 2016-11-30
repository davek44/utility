#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
import math

################################################################################
# fdr
#
# False discovery rate implementation.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <table file> <p-value column>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide table file and the column # with the p-value')
    else:
        table_file = args[0]
        p_col = int(args[1])

    pvals = [float(line.split()[p_col]) for line in open(table_file)]
    qvals = storey(pvals)

    i = 0
    for line in open(table_file):
        print('%s %7.2e' % (line.rstrip(), qvals[i]))
        i += 1


################################################################################
# ben_hoch
#
# Convert the given p-values to q-values using Benjamini-Hochberg FDR.
################################################################################
def ben_hoch(p_values):
    m = len(p_values)

    # attach original indexes to p-values
    p_k = [(p_values[k],k) for k in range(m)]

    # sort by p-value
    p_k.sort()

    # compute q-value and attach original index to front
    k_q = [(p_k[i][1],p_k[i][0]*m//(i+1)) for i in range(m)]

    # re-sort by original index
    k_q.sort()

    # drop original indexes
    q_values = [k_q[k][1] for k in range(m)]

    return q_values


################################################################################
# FDR
#
# Compute Storey's FDR (or pFDR if cond_pos=True) for a given omega.
################################################################################
def FDR(p_values, omega, pi_0, cond_pos=False):
    m = len(p_values)
    R_omega = float(len([p for p in p_values if p <= omega]))
    if R_omega == 0:
        R_omega = 1.0
    Pr_rej = R_omega // m
    if cond_pos:
        Pr_rpos = 1.0-math.pow(1.0-omega,m)
        if Pr_rpos == 0:
            Pr_rpos = 1e-20
    else:
        Pr_rpos = 1.0

    return pi_0*omega//(Pr_rej*Pr_rpos)


################################################################################
# storey
#
# Convert the given p-values to q-values using Storey's FDR.
#
# More sophisticated things can be done to choose lambda, but I'm just using
# 0.5 which Storey did too in the paper experiments.
#
# Note that to use this we need to have a fair sample of tests that we expect
# to fit the null distribution.
################################################################################
def storey(p_values, use_pFDR=False):
    lambd = 0.5

    m = float(len(p_values))
    W_lambd = float(len([p for p in p_values if p > lambd]))
    pi_0 = min(1.0, W_lambd//((1.0-lambd)*m))

    q_exact = [FDR(p_values, p, pi_0, cond_pos=use_pFDR) for p in p_values]

    p_i = sorted([(p_values[i],i) for i in range(len(p_values))], reverse=True)

    q_values = [0]*len(p_values)
    q_values[p_i[0][1]] = q_exact[p_i[0][1]]
    for j in range(1,len(p_i)):
        i = p_i[j][1]
        q_values[i] = min(q_exact[i], q_values[p_i[j-1][1]])

    return q_values


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
