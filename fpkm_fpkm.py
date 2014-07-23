#!/usr/bin/env python
from optparse import OptionParser
import os, math, sys
from scipy.stats import spearmanr
import gff, ggplot, cufflinks

################################################################################
# fpkm_fpkm.py
#
# Compare two cufflinks runs.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fpkm1_file> <fpkm2_file>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='gtf')
    parser.add_option('-o', dest='out_dir', default='.')
    parser.add_option('-p', dest='pseudocount', default=0.125, type='float')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide two diff files')
    else:
        fpkm1_file = args[0]
        fpkm2_file = args[1]

    cuff1 = cufflinks.fpkm_tracking(fpkm1_file)
    cuff2 = cufflinks.fpkm_tracking(fpkm2_file)

    gtf_genes = set()
    if options.gtf:
        gtf_genes = gff.gtf_gene_set(options.gtf)

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    for sample in cuff1.experiments:
        # scatter plot fpkm
        df = {'fpkm1':[], 'fpkm2':[]}
        for i in range(len(cuff1.genes)):
            if len(gtf_genes) == 0 or cuff1.genes[i] in gtf_genes:
                fpkm1 = cuff1.gene_expr_exp(i, sample)
                fpkm2 = cuff2.gene_expr_exp(i, sample)

                if not math.isnan(fpkm1) and not math.isnan(fpkm2):
                    df['fpkm1'].append(math.log(options.pseudocount+fpkm1,2))
                    df['fpkm2'].append(math.log(options.pseudocount+fpkm2,2))

        r_script = '%s/fpkm_fpkm_scatter.r' % os.environ['RDIR']
        out_pdf = '%s/%s_scatter.pdf' % (options.out_dir, sample)
        ggplot.plot(r_script, df, [out_pdf])

        # compute correlation
        cor, p = spearmanr(df['fpkm1'], df['fpkm2'])

        report_out = open('%s/%s_report.txt' % (options.out_dir,sample), 'w')
        print >> report_out, 'Spearman correlation: %f (%e)' % (cor, p)
        report_out.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
