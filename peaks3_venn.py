#!/usr/bin/env python
from optparse import OptionParser
import math, os, stats, subprocess, tempfile

import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import seaborn as sns

################################################################################
# peaks3_venn.py
#
# Make a venn diagram comparing three sets of peak calls.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <peaks1_bed> <peaks2_bed> <peaks3_bed> <out_pdf>'
    parser = OptionParser(usage)
    parser.add_option('--l1', dest='label1', default='peaks1', help='Label for peak set 1')
    parser.add_option('--l2', dest='label2', default='peaks2', help='Label for peak set 2')
    parser.add_option('--l3', dest='label3', default='peaks3', help='Label for peak set 3')
    (options,args) = parser.parse_args()

    if len(args) != 4:
        parser.error('Must provide three peaks BED files and output PDF')
    else:
        peak_beds = args[:3]
        out_pdf = args[3]

    merge_fd, merge_bed = tempfile.mkstemp()

    # merge peaks
    cmd = 'cat %s %s %s | awk \'{OFS="\t"} {print $1, $2, $3}\' | bedtools sort -i stdin | bedtools merge -i stdin > %s' % (peak_beds[0], peak_beds[1], peak_beds[2], merge_bed)
    subprocess.call(cmd, shell=True)

    # annotate merged peaks with each individual set
    num_peaks = count_peaks(merge_bed)
    peak_overlaps = [set(), set(), set()]

    for bi in range(3):
        cmd = 'bedtools intersect -c -a %s -b %s' % (merge_bed, peak_beds[bi])
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        pi = 0
        for line in p.stdout:
            a = line.split()
            if int(a[-1]) > 0:
                peak_overlaps[bi].add(pi)
            pi += 1

    # plot
    plt.figure()
    venn_diag = venn3(peak_overlaps, set_labels=[options.label1, options.label2, options.label3]) # , set_colors=['#e41a1c', '#A1A838', ''])
    plt.savefig(out_pdf)
    plt.close()

    # clean up
    os.close(merge_fd)
    os.remove(merge_bed)


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
