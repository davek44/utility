#!/usr/bin/env python
from optparse import OptionParser
import sys
import gff

################################################################################
# gtf_cut.py
#
# Cut a gtf key:value pair out of a gtf file.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] -k <key> <gtf file>'
    parser = OptionParser(usage)
    parser.add_option('-k', dest='key', help='Key to extract')
    (options,args) = parser.parse_args()

    if len(args) == 1:
        if args[0] == '-':
            gtf_open = sys.stdin
        else:
            gtf_open = open(args[0])
    else:
        parser.error(usage)

    if not options.key:
        parser.error('Must provide key')

    for line in gtf_open:
        a = line.split('\t')
        kv = gff.gtf_kv(a[8])
        print kv.get(options.key,'')

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
