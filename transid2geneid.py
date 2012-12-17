#!/usr/bin/env python
from optparse import OptionParser
import gff

################################################################################
# transid2geneid.py
#
# Given a transcript id, produce a gene id to punch into the browser
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <trans id>'
    parser = OptionParser(usage)
    parser.add_option('-l', dest='lnc_file', default='/Users/dk/research/common/data/lncrna/lnc_catalog.gtf', help='lncRNA catalog file [Default: %default]')
    (options,args) = parser.parse_args()
    
    if len(args) != 1:
        parser.error('Must provide transcript id')
    else:
        trans_id = args[0]

    for line in open(options.lnc_file):
        a = line.split('\t')
        kv = gff.gtf_kv(a[8])
        if kv['transcript_id'] == trans_id:
            print kv['gene_id']
            break


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
