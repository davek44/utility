#!/usr/bin/env python
from optparse import OptionParser
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.lib.ggplot2 as ggplot2
import pdb
import pysam

grdevices = importr('grDevices')

################################################################################
# bam_len_hist.py
#
# Plot a histogram of the length of alignments in a BAM file.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide BAM file')
    else:
        bam_file = args[0]

    align_lengths = {}
    for aligned_read in pysam.Samfile(bam_file, 'rb'):
        align_lengths[aligned_read.qlen] = align_lengths.get(aligned_read.qlen,0) + 1

    min_len = min(align_lengths.keys())
    max_len = max(align_lengths.keys())

    # construct data frame
    len_r = ro.IntVector(range(min_len,max_len+1))
    counts_r = ro.IntVector([align_lengths.get(l,0) for l in range(min_len,max_len+1)])
    
    df = ro.DataFrame({'length':len_r, 'counts':counts_r})

    # construct full plot
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='length', y='counts') + \
        ggplot2.geom_bar(stat='identity') + \
        ggplot2.scale_x_continuous('Alignment length') + \
        ggplot2.scale_y_continuous('')

    # plot to file
    grdevices.pdf(file='align_lengths.pdf')
    gp.plot()
    grdevices.dev_off()

    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
