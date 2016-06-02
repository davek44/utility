#!/usr/bin/env python
from optparse import OptionParser
import numpy as np
import math

################################################################################
# ripseq.py
#
# Methods to aid RIP-Seq analysis.
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
# hash_input_fpkm
################################################################################
def hash_input_fpkm(diff_file, just_ok=False, log_pseudo=None, by_rbp=False):
    input_fpkm_raw = {}

    diff_in = open(diff_file)
    diff_in.readline()
    for line in diff_in:
        a = line.split('\t')

        gene_id = a[0]
        sample1 = a[4]
        sample2 = a[5]
        status = a[6]
        fpkm1 = float(a[7])
        fpkm2 = float(a[8])

        if sample1 == 'input':
            rbp = sample2
            ifpkm = fpkm1
        elif sample2 == 'input':
            rbp = sample1
            ifpkm = fpkm2
        else:
            print >> sys.stderr, 'Cannot find input: %s' % line,

        if not just_ok or status == 'OK':
            if log_pseudo == None:
                fpkm_final = ifpkm
            else:
                fpkm_final = np.log2(ifpkm+log_pseudo)

            if by_rbp:
                input_fpkm_raw.setdefault(rbp,{})[gene_id] = fpkm_final
            else:
                #input_fpkm[gene_id] = fpkm_final
                input_fpkm_raw.setdefault(gene_id,[]).append(fpkm_final)

    diff_in.close()

    # average multiples
    if by_rbp:
        input_fpkm = input_fpkm_raw

    else:
        input_fpkm = {}
        for gene_id in input_fpkm_raw:
            input_fpkm[gene_id] = np.mean(input_fpkm_raw[gene_id])

    return input_fpkm


################################################################################
# hash_rip_fold
#
# Input:
#  diff_file:    RIP gene_exp.diff file, comparing RBPs and input.
#  min_fpkm:     Minimum FPKM to consider a gene.
#  pseudocount:  Pseudocount to allow safer fold changes for zero FPKM genes.
#  max_fold:     Max to cap fold change at.
#  one_rbp:      diff_file contains only a single RBP, so return a hash of
#                   gene_id to fold changes instead.
#
# Output:
#  rbp_fold:     Dict mapping RBPs to dicts mapping gene_id to fold changes
################################################################################
def hash_rip_fold(diff_file, min_fpkm=0, pseudocount=0.125, max_fold=None, one_rbp=False):
    rbp_fold = {}
    rbp_bound = {}

    # read rip diff
    diff_in = open(diff_file)
    diff_in.readline()
    for line in diff_in:
        a = line.split('\t')

        gene_id = a[0]
        sample1 = a[4]
        sample2 = a[5]
        fpkm1 = float(a[7])
        fpkm2 = float(a[8])

        rbp = None
        if sample1 == 'input':
            rbp = sample2
        elif sample2 == 'input':
            rbp = sample1
            fpkm1, fpkm2 = fpkm2, fpkm1

        # if fpkm1 > min_fpkm or fpkm2 > min_fpkm:
        if fpkm1 > min_fpkm and fpkm2 > min_fpkm:
            rip_fold = np.log2(fpkm2+pseudocount) - np.log2(fpkm1+pseudocount)

            if max_fold != None:
                rip_fold = min(rip_fold, max_fold)
                rip_fold = max(rip_fold, -max_fold)

            if one_rbp:
                rbp_fold[gene_id] = rip_fold
            else:
                rbp_fold.setdefault(rbp,{})[gene_id] = rip_fold

    diff_in.close()

    return rbp_fold


################################################################################
# hash_rip
#
# Input:
#  diff_file:      RIP gene_exp.diff file, comparing RBPs and input.
#  just_ok:        Return the RIP stat for only status OK genes.
#  min_fpkm:       Min FPKM required in either sample to consider a gene.
#  max_stat:       Max to cap the RIP statistic at.
#  use_fold:       Use fold change rather than test_stat as the RIP stat.
#  one_rbp:        There is only one RBP, so hash genes directly.
#
# Output:
#  rbp_gene_tstat: Dict mapping RBPs to dicts mapping gene_id to test_stat
#  rbp_bound:      Dict mapping RBPs to sets of bound gene_ids.
################################################################################
def hash_rip(diff_file, just_ok=False, use_fold=False, max_stat=None, min_fpkm=None, pseudocount=0.125, one_rbp=False):
    rbp_stat = {}
    rbp_bound = {}

    # read rip diff
    diff_in = open(diff_file)
    diff_in.readline()
    for line in diff_in:
        a = line.split('\t')

        gene_id = a[0]
        sample1 = a[4]
        sample2 = a[5]
        status = a[6]
        fpkm1 = float(a[7])
        fpkm2 = float(a[8])
        fold_change = float(a[9])
        tstat = float(a[10])
        qval = float(a[11])
        sig = a[-1].rstrip()

        rbp = None
        if sample1 == 'input':
            rbp = sample2
        elif sample2 == 'input':
            rbp = sample1
            fpkm1, fpkm2 = fpkm2, fpkm1
            fold_change *= -1
            tstat *= -1

        if (not just_ok or status == 'OK') and not math.isnan(tstat):
            if min_fpkm == None or fpkm1 > min_fpkm or fpkm2 > min_fpkm:
                if use_fold:
                    # rip_stat = fold_change
                    rip_stat = np.log2(fpkm2+pseudocount) - np.log2(fpkm1+pseudocount)
                else:
                    rip_stat = tstat

                if max_stat != None:
                    rip_stat = min(rip_stat, max_stat)
                    rip_stat = max(rip_stat, -max_stat)

                if one_rbp:
                    rbp_stat[gene_id] = rip_stat
                else:
                    rbp_stat.setdefault(rbp,{})[gene_id] = rip_stat

                if sig == 'yes':
                    if tstat > 0:
                        if one_rbp:
                            rbp_bound[gene_id] = True
                        else:
                            rbp_bound.setdefault(rbp,{})[gene_id] = True
                    else:
                        if one_rbp:
                            rbp_bound[gene_id] = False
                        else:
                            rbp_bound.setdefault(rbp,{})[gene_id] = False

    diff_in.close()

    # add rbps w/ no sig genes
    if not one_rbp:
        for rbp in rbp_stat:
            if not rbp in rbp_bound:
                rbp_bound[rbp] = set()

    return rbp_stat, rbp_bound


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
