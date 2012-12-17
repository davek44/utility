#!/usr/bin/env python
from optparse import OptionParser

################################################################################
# fastq_quality_change.py
#
# Change the quality value ascii index for a fastq file.
#
# Author: David Kelley dakelley@umiacs.umd.edu
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fastq file>'
    parser = OptionParser(usage)
    parser.add_option('-a', dest='after', type='int', help='Desired fastq quality ascii index')
    parser.add_option('-b', dest='before', type='int', help='Current fastq quality ascii index')
    (options,args) = parser.parse_args()

    if options.after == None or options.before == None:
        parser.error('Must provide before and after ascii indexes')
    if len(args) != 1:
        parser.error('Must provide fastq file')
    else:
        fastq_file = args[0]

    fq_in = open(fastq_file)
    header = fq_in.readline()
    while header:        
        seq = fq_in.readline()
        mid = fq_in.readline()
        qual = fq_in.readline()

        print header,
        print seq,
        print mid,
        print ''.join([chr(ord(q)-options.before+options.after) for q in qual])

        header = fq_in.readline()
    fq_in.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
