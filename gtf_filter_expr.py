#!/usr/bin/env python
from optparse import OptionParser
import cufflinks, gff

################################################################################
# gtf_filter_expr.py
#
# Filter a gtf file to only leave genes that are expressed over a specified
# threshold in a specified tissue/cell type.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf file> <cell type>'
    parser = OptionParser(usage)
    parser.add_option('-t', dest='expr_t', type='float', default=.1, help='Minimum allowed fpkm value')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide gtf file and cell type')
    else:
        gtf_file = args[0]
        cell_type = args[1]

    # get expression data
    cuff = cufflinks.fpkm_tracking()

    # find cell type experiment index
    cell_indexes = [i for i in range(len(cuff.experiments)) if cuff.experiments[i]==cell_type]
    if len(cell_indexes) == 0:
        parser.error('Cell type %s does not match any quantified experiments' % cell_type)
    else:
        cell_i = cell_indexes[0]

    # parser gtf file
    for line in open(gtf_file):
        a = line.split('\t')
        gene_id = gff.gtf_kv(a[8])['gene_id']
        expr_vec = cuff.gene_expr(gene_id)
        if expr_vec[cell_i] > options.expr_t:
            print line,
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
