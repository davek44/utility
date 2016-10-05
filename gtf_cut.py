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
    parser.add_option('-l', dest='line_too', action='store_true', default=False, help='Print the line too [Default: %default]')
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
    else:
        keys = options.key.split(',')

    for line in gtf_open:
        if not line.startswith('##'):
            a = line.split('\t')
            kv = gff.gtf_kv(a[8])

            if options.line_too:
                key_str = '\t'.join([kv.get(key,'-') for key in keys])
                print('%s\t%s' % (key_str,line))
            else:
                print('\t'.join([kv.get(key,'-') for key in keys]))


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
