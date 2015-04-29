#!/usr/bin/env python
from optparse import OptionParser
import math, os

################################################################################
# gsea_rnk.py
#
# Output a set of .rnk files for GSEA from a cuffdiff .diff file.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <diff>'    
    parser = OptionParser(usage)
    parser.add_option('-m', dest='min_fpkm', type='float')
    parser.add_option('-o', dest='out_dir', default='.')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide .diff')
    else:
        diff_file = args[0]

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    comparison_out = {}

    diff_in = open(diff_file)
    diff_in.readline()
    for line in diff_in:
        a = line.split('\t')

        gene_id = a[0]
        gene_name = a[2]
        sample1 = a[4].replace('-','_')  # cmd line gsea cannot handle hyphens
        sample2 = a[5].replace('-','_')
        status = a[6]
        fpkm1 = float(a[7])
        fpkm2 = float(a[8])
        fold_change = float(a[9])
        tstat = float(a[10])
        qval = float(a[11])
        sig = a[-1].rstrip()

        if status == 'OK' and not math.isnan(tstat):
            if options.min_fpkm == None or fpkm1 > options.min_fpkm or fpkm2 > options.min_fpkm:
                if not (sample1,sample2) in comparison_out:
                    comparison_out[(sample1,sample2)] = open('%s/%s_%s.rnk' % (options.out_dir, sample1, sample2), 'w')

                print >> comparison_out[(sample1,sample2)], '%s\t%f' % (gene_name, fold_change)

    diff_in.close()

    for ckey in comparison_out:
        comparison_out[ckey].close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
