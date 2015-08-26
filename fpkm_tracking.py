#!/usr/bin/env python
from optparse import OptionParser
import pandas as pd

################################################################################
# fpkm_tracking
#
# Print a table of FPKM abundance estimates with one gene/sample per row.
#
# Using Pandas for this is stupid because we have to read the whole thing
# into memory, and when you iterate over rows, it has to create a Series object.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fpkm_tracking>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='gene_id', help='This gene only')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide .fpkm_tracking file')
    else:
        fpkm_tracking_file = args[0]
    
    cuff = pd.read_csv(fpkm_tracking_file, sep='\t')

    fpkm_indexes = [i for i in range(cuff.shape[1]) if cuff.columns[i][-5:] == '_FPKM']

    for gene_i, gene_series in cuff.iterrows():
        gene_id = gene_series['gene_id']
        if options.gene_id == None or gene_id == options.gene_id:
            for i in fpkm_indexes:
                sample = cuff.columns[i][:-5]
                fpkm = str(gene_series[i])
                status = gene_series[i+3]

                if status == 'OK':
                    cols = [gene_series['tracking_id'], sample, fpkm]
                    print '\t'.join(cols)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
