#!/usr/bin/env python
from optparse import OptionParser
import math
import numpy as np

################################################################################
# cuffdiff.py
#
# Methods to aid analysis of cuffdiff output
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
# hash_stat
# 
# Input:
#  diff_file:     RIP *_exp.diff file.
#  stat:          test_stat or fold
#  max_stat:      Maximum abs value allowed for the diff stat.
#  min_fpkm:      Minimum FPKM to consider a gene.
#  pseudocount:   Pseudocount to compute my own log fold changes.
#  sample_first:  Sample name to force to come first (e.g input/control/etc)
#  gene_set:      Only return genes in this set.
#
# Output:
#  gene_diff:     Dict mapping sample pairs to dicts mapping gene_id to diff stat
################################################################################
def hash_stat(diff_file, stat='fold', max_stat=None, min_fpkm=None, pseudocount=0.125, sample_first=None, gene_set=None):
    gene_diff = {}

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
        # fold_change = float(a[9])
        test_stat = float(a[10])
        qval = float(a[11])
        sig = a[-1].rstrip()

        if sample2 == sample_first:
            sample2, sample1 = sample1, sample2
            fpkm2, fpkm1 = fpkm1, fpkm2
            # fold_change *= -1
            test_stat *= -1

        if gene_set is None or gene_id in gene_set:
            if min_fpkm is None or fpkm1 > min_fpkm or fpkm2 > min_fpkm:
                if stat in ['fold','fold_change']:
                    diff_stat = np.log2(fpkm2+pseudocount) - np.log2(fpkm1+pseudocount)
                elif stat in ['test_stat','tstat']:
                    if status == 'OK' and not math.isnan(test_stat):
                        diff_stat = test_stat
                    else:
                        diff_stat = None
                else:
                    print >> sys.stderr, 'Unknown stat requested: %s' % stat
                    exit(1)

                if diff_stat is not None:
                    if max_stat is not None:
                        diff_stat = min(diff_stat, abs(max_stat))
                        diff_stat = max(diff_stat, -abs(max_stat))

                    gene_diff.setdefault((sample1,sample2),{})[gene_id] = diff_stat

    diff_in.close()

    return gene_diff


################################################################################
# hash_sig
# 
# Input:
#  diff_file:     RIP *_exp.diff file.
#  sample_first:  Sample name to force to come first (e.g input/control/etc)
#  gene_set:      Only return genes in this set.
#
# Output:
#  gene_diff:     Dict mapping sample pairs to dicts mapping gene_id to diff stat
################################################################################
def hash_sig(diff_file, sample_first=None, gene_set=None):
    gene_sig = {}

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
        test_stat = float(a[10])
        qval = float(a[11])
        sig = a[-1].rstrip()

        if sample2 == sample_first:
            sample2, sample1 = sample1, sample2
            fpkm2, fpkm1 = fpkm1, fpkm2
            fold_change *= -1
            test_stat *= -1

        if gene_set is None or gene_id in gene_set:
            if (sample1,sample2) not in gene_sig:
                gene_sig[(sample1,sample2)] = {}

            if sig == 'yes':
                if test_stat > 0:
                    sig_val = 1
                else:
                    sig_val = -1
                gene_sig[(sample1,sample2)][gene_id] = sig_val

    diff_in.close()

    return gene_sig


################################################################################
# hash_stat_one
# 
# Input:
#  diff_file:     RIP *_exp.diff file.
#  stat:          test_stat or fold
#  max_stat:      Maximum abs value allowed for the diff stat.
#  min_fpkm:      Minimum FPKM to consider a gene.
#  sample_first:  Sample name to force to come first (e.g input/control/etc)
#  gene_set:      Only return genes in this set.
#
# Output:
#  gene_diff:     Dict mapping gene_id to diff stat
################################################################################
def hash_stat_one(diff_file, stat='fold', max_stat=None, min_fpkm=None, pseudocount=0.125, sample_first=None):
    gene_diff = hash_stat(diff_file, stat, max_stat, min_fpkm, pseudocount, sample_first, gtf_genes)
    if len(gene_diff.keys()) > 1:
        print >> sys.stderr, 'More than one pair of samples found in %s' % diff_file
        exit(1)
    else:
        samples = gene_diff.keys()[0]
        return gene_diff[samples]


################################################################################
# hash_sig_one
# 
# Input:
#  diff_file:     RIP *_exp.diff file.
#  sample_first:  Sample name to force to come first (e.g input/control/etc)
#  gene_set:      Only return genes in this set.
#
# Output:
#  gene_diff:     Dict mapping gene_id to diff stat
################################################################################
def hash_sig_one(diff_file, sample_first=None, gtf_genes=None):
    gene_sig = hash_sig(diff_file, sample_first, gtf_genes)
    if len(gene_sig.keys()) > 1:
        print >> sys.stderr, 'More than one pair of samples found in %s' % diff_file
        exit(1)
    else:
        samples = gene_sig.keys()[0]
        return gene_sig[samples]


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
