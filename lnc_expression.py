#!/usr/bin/env python
from optparse import OptionParser
import cufflinks, gff
import os

################################################################################
# lnc_expession.py
#
# Print a summary of the lncrna gene's expression.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gene/transcript id>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='cuff_dir', default='%s/research/common/data/lncrna'%os.environ['HOME'], help='Cufflinks output directory with .fpkm_tracking files [Default: %default]')
    parser.add_option('-l', dest='lnc_gtf', default='%s/research/common/data/lncrna/lnc_catalog.gtf'%os.environ['HOME'], help='lncRNA catalog gtf file [Default: %default]')
    parser.add_option('-t', dest='transcript_expr', default=False, action='store_true', help='Return transcript expression rather than gene [Default: %default]')
    (options,args) = parser.parse_args()

    if options.transcript_expr:
        cuff = cufflinks.fpkm_tracking('%s/isoforms.fpkm_tracking' % options.cuff_dir)

        if args[0].find('XLOC') != -1:
            trans_ids = set()
            for line in open(options.lnc_gtf):
                a = line.split('\t')
                kv = gff.gtf_kv(a[8])
                if kv['gene_id'] == args[0]:
                    trans_ids.add(kv['transcript_id'])
        else:
            trans_ids = [args[0]]

        for trans_id in trans_ids:
            print '%s:' % trans_id
            cuff.gene_expr_print(trans_id)

    else:
        cuff = cufflinks.fpkm_tracking('%s/genes.fpkm_tracking' % options.cuff_dir)

        if args[0].find('XLOC') != -1:
            gene_id = args[0]
        else:
            t2g = gff.t2g(options.lnc_gtf)
            gene_id = t2g[args[0]]

        cuff.gene_expr_print(gene_id)
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
