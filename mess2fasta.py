#!/usr/bin/env python
from optparse import OptionParser

################################################################################
# mess2fa.py
#
# Convert a mess of nt's and numbers and spaces into a neat fasta file
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <mess file>'
    parser = OptionParser(usage)
    parser.add_option('-d','--header', dest='header', default='mess', help='Fasta header [Default: %default]')
    parser.add_option('-u', '--upper', dest='upper', action='store_true', default=False, help='Uppercase all nucleotides [Default: %default]')
    (options,args) = parser.parse_args()

    allowed_nts = set(['A','C','G','T','N','a','c','g','t','n'])

    seq = ''
    for line in open(args[0]):
        seq += ''.join([nt for nt in line if nt in allowed_nts])

    if options.upper:
        seq = ''.join([nt.upper() for nt in seq])

    print '>%s\n%s' % (options.header,seq)
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
