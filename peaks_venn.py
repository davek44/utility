#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess
import math, os, stats, subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

################################################################################
# peaks_venn.py
#
# Make a venn diagram comparing two sets of peak calls.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <peaks1_bed> <peaks2_bed> <out_pdf>'
    parser = OptionParser(usage)
    parser.add_option('--l1', dest='label1', default='peaks1', help='Label for peak set 1')
    parser.add_option('--l2', dest='label2', default='peaks2', help='Label for peak set 2')
    (options,args) = parser.parse_args()

    if len(args) != 3:
        parser.error('Must provide two peaks BED files and output PDF')
    else:
        peaks1_bed = args[0]
        peaks2_bed = args[1]
        out_pdf = args[2]

    # count individual
    peaks1_count = count_peaks(peaks1_bed)
    peaks2_count = count_peaks(peaks2_bed)

    # count overlap
    copeaks_count = 0
    p = subprocess.Popen('intersectBed -u -a %s -b %s' % (peaks1_bed, peaks2_bed), stdout=subprocess.PIPE, shell=True)
    for line in p.stdout:
        copeaks_count += 1
    p.communicate()

    plt.figure()
    venn_diag = venn2(subsets=(peaks1_count-copeaks_count, peaks2_count-copeaks_count, copeaks_count), set_labels=[options.label1, options.label2], set_colors=['#e41a1c', '#A1A838'])
    plt.savefig(out_pdf)
    plt.close()


def count_peaks(bed_file):
    peak_counts = 0
    for line in open(bed_file):
        peak_counts += 1
    return peak_counts

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
