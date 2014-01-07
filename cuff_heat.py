#!/usr/bin/env python
from optparse import OptionParser
import math, os, pdb, random, sys
import cufflinks, gff, ggplot

################################################################################
# cuff_heat.py
#
################################################################################

log_min = .25

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fpkm_tracking>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='gtf', help='GTF file of genes to display')
    parser.add_option('-m', dest='min_fpkm', default=.125, help='Minimum FPKM (for logs) [Default: %default]')
    parser.add_option('-o', dest='out_pdf', default='cuff_heat.pdf', help='Output PDF [Default: %default]')
    parser.add_option('-s', dest='sample', default=50, help='Sample genes rather than use all [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide fpkm_tracking')
    else:
        fpkm_tracking = args[0]

    # load expression data
    cuff = cufflinks.fpkm_tracking(fpkm_file=fpkm_tracking)

    # determine genes
    all_genes = cuff.genes
    if options.gtf:
        all_genes = set()
        for line in open(options.gtf):
            a = line.split('\t')
            all_genes.add(gff.gtf_kv(a[8])['gene_id'])

    # sample genes to display
    if len(all_genes) <= options.sample:
        display_genes = all_genes
    else:
        display_genes = random.sample(all_genes, options.sample)

    # build data frame
    df = {'Gene':[], 'FPKM':[], 'Sample':[]}

    for gene_id in display_genes:
        ge = cuff.gene_expr(gene_id)
        if not math.isnan(ge[0]):
            for i in range(len(cuff.experiments)):
                df['Gene'].append(gene_id)
                df['Sample'].append(cuff.experiments[i])
                df['FPKM'].append(math.log(ge[i]+options.min_fpkm,2))

    # plot
    ggplot.plot('%s/cuff_heat.r' % os.environ['RDIR'], df, [options.out_pdf], debug=True)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
