#!/usr/bin/env python
from optparse import OptionParser

################################################################################
# cuff_fails.py
#
# Print a table of the number of genes for which cufflinks failed.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fpkm_tracking>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide FPKM tracking file')
    else:
        fpkm_file = args[0]

    # get headers
    fpkm_in = open(fpkm_file)
    headers = fpkm_in.readline().split()

    # initialize fail counts
    fails = {}

    for line in fpkm_in:
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        gene_id = a[0]

        for i in range(len(a)):
            if headers[i][-7:] == '_status':
                if a[i] != 'OK':
                    sample = headers[i][:-7]
                    fails[sample] = fails.get(sample,0) + 1
        
    fpkm_in.close()

    for sample in fails:
        cols = (sample, fails[sample])
        print '%-18s  %5d' % cols


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
