#!/usr/bin/env python
from optparse import OptionParser
import numpy as np
import math, os, pdb, random, sys
import ggplot

################################################################################
# cuff_bar.py
#
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fpkm_tracking> <gene1> <gene2> ...>'
    parser = OptionParser(usage)
    parser.add_option('-l', dest='log', default=False, help='log2 FPKM')
    parser.add_option('-n', dest='names', default=None, help='Sample names, comma-separated')
    parser.add_option('-p', dest='pseudocount', default=.125, type='float', help='Pseudocount for log FPKM [Default: %default]')
    parser.add_option('-o', dest='out_dir', default='cuff_bars', help='Output directory [Default: %default]')
    parser.add_option('-s', dest='samples', default=None, help='Samples to plot, comma-separated')
    parser.add_option('-y', dest='yaxis_match', default=False, action='store_true', help='Match the y-axis of all plots [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) < 2:
        parser.error('Must provide fpkm_tracking and genes.')
    else:
        fpkm_tracking = args[0]
        genes = args[1:]

    gene_sample_fpkm = read_fpkms(fpkm_tracking, genes, options.log, options.pseudocount)

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    if options.samples:
        samples = options.samples.split(',')
    else:
        samples = sorted(gene_sample_fpkm[genes[0]].keys())

    if options.names:
        names = options.names.split(',')
    else:
        names = samples

    ymin = 0
    if options.log:
        ymin = np.log2(options.pseudocount)

    if options.yaxis_match:
        ymax = max([gene_sample_fpkm[gene_name][sample][2] for sample in samples for gene_name in genes])
    else:
        ymax = None

    for gene_name in genes:
        df = {}
        df['Sample'] = names
        df['FPKM'] = [gene_sample_fpkm[gene_name][sample][0] for sample in samples]
        df['conf_lo'] = [gene_sample_fpkm[gene_name][sample][1] for sample in samples]
        df['conf_hi'] = [gene_sample_fpkm[gene_name][sample][2] for sample in samples]

        out_pdf = '%s/%s.pdf' % (options.out_dir, gene_name)
        ggplot.plot('%s/cuff_bar.r'%os.environ['RDIR'], df, [ymin, ymax, out_pdf])


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
# read_fpkms
################################################################################
def read_fpkms(fpkm_file, genes, log_fpkm, pseudo):
    gene_sample_fpkm = {}

    fpkm_in = open(fpkm_file)
    
    # parse headers into sample names
    headers = fpkm_in.readline().split()
    samples = [h[:-5] for h in headers if h == 'FPKM' or h[-5:] == '_FPKM']

    for line in fpkm_in:
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        gene_name = a[4]
        if gene_name in genes:
            gene_sample_fpkm[gene_name] = {}

            for i in range(len(headers)):
                if headers[i] == 'FPKM' or headers[i][-5:] == '_FPKM':
                    sample = headers[i][:-5]
                    status = a[i+3]

                    if status in ['FAIL','HIDATA']:
                        if log_fpkm:
                            gene_sample_fpkm[gene_name][sample] = tuple(3*[np.log2(pseudo)])
                        else:
                            gene_sample_fpkm[gene_name][sample] = (0, 0, 0)
                    else:
                        if log_fpkm:
                            gene_sample_fpkm[gene_name][sample] = (np.log2(pseudo+float(a[i])), np.log2(pseudo+float(a[i+1])), np.log2(pseudo+float(a[i+2])))
                        else:
                            gene_sample_fpkm[gene_name][sample] = (float(a[i]), float(a[i+1]), float(a[i+2]))

    fpkm_in.close()

    return gene_sample_fpkm

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
