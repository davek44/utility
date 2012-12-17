#!/usr/bin/env python
from optparse import OptionParser
import sys
import gff

################################################################################
# gtf_filter_csf.py
#
# Filter the lnc catalog gtf file by CSF value.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf file>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='greater', action='store_true', default=False, help='Keep genes w/ CSF value greater than the one given [Default: %default]')
    parser.add_option('-l', dest='less', action='store_true', default=True, help='Keep genes w/ CSF value less than the one given [Default: %default]')
    parser.add_option('-t', dest='csf_t', type='float', default=100.0, help='CSF threshold [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) == 1:
        gtf_open = open(args[0])
    else:
        gtf_open = sys.stdin

    line = gtf_open.readline()
    while line:
        a = line.split('\t')
        csf = float(gff.gtf_kv(a[8])['csf'])
        if (options.less and csf <= options.csf_t) or (options.greater and csf >= options.csf_t):
            print line,
        line = gtf_open.readline()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
