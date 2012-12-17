#!/usr/bin/env python
from optparse import OptionParser
import math, os, pdb, random
import cufflinks, gff

from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.lib.ggplot2 as ggplot2
grdevices = importr('grDevices')

################################################################################
# fpkm_hist.py
#
# Plot a histogram of the max log2 FPKM values for the genes in a gtf file.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf file> <fpkm tracking>'
    parser = OptionParser(usage)
    #parser.add_option('-m', dest='fpkm_min', type='float', default=0.25, help='Minimum FPKM [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error(usage)
    else:
        gtf_file = args[0]
        fpkm_tracking_file = args[1]

    # get genes
    genes = set()
    for line in open(gtf_file):
        a = line.split('\t')
        genes.add(gff.gtf_kv(a[8])['gene_id'])

    # get expression
    cuff = cufflinks.fpkm_tracking(fpkm_tracking_file)
    log_fpkms = []
    for gene_id in genes:
        max_fpkm = max(cuff.gene_expr(gene_id))
        if max_fpkm > 0:
            log_fpkms.append(math.log(max_fpkm,2))

    # construct R data objects
    fpkms_r = ro.FloatVector(log_fpkms)
    df = ro.DataFrame({'fpkm':fpkms_r})
    
    # construct plot
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='fpkm') + \
        ggplot2.geom_histogram(binwidth=0.2)
    
    # save to file
    gtf_pre = os.path.splitext(gtf_file)[0]
    grdevices.pdf(file='%s_fpkmhist.pdf' % gtf_pre)
    gp.plot()
    grdevices.dev_off()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
