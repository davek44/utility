#!/usr/bin/env python
from optparse import OptionParser
import pysam

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import seaborn as sns
sns.set_style('ticks')

################################################################################
# plot_fragment_lengths.py
#
# Plot the distribution of fragment lengths
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam> <out_pdf>'
    parser = OptionParser(usage)
    parser.add_option('-m', dest='max_length', type='int')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide BAM file and output PDF')
    else:
        bam_file = args[0]
        out_pdf = args[1]

    tlens = {}
    for alignment in pysam.Samfile(bam_file):
        tl = abs(alignment.template_length)
        tlens[tl] = tlens.get(tl,0) + 1

    # not sure what 0 means
    tlens[0] = 0

    if options.max_length is None:
        num_fragments = sum([tlens.get(i,0) for i in range(10000)])
        max_length_fragments = 0.99*num_fragments

        length = 1
        length_fragments = tlens.get(length,0)
        while length_fragments < max_length_fragments and length < 1000:
            print length, length_fragments, max_length_fragments
            length += 1
            length_fragments += tlens.get(length,0)

        options.max_length = length

    #for tl in range(max(tlens.keys())+1):
    #    print '%4d  %d' % (tl,tlens.get(tl,0))

    length_counts = [tlens.get(length,0) for length in range(options.max_length)]

    plt.figure()
    plt.plot(length_counts)
    plt.xlabel('Fragment length')
    plt.xlim(0,options.max_length+1)
    sns.despine()
    plt.savefig(out_pdf)
    plt.close()
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
