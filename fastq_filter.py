#!/usr/bin/env python
from optparse import OptionParser
import bz2
import gzip

'''
fastq_filter.py

Filter a FASTQ file for various properties, like read length.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fastq_file>'
    parser = OptionParser(usage)
    parser.add_option('-l', dest='length_min', default=None, type='int', help='Minimum read length')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide FASTQ file')
    else:
        fastq_file = args[0]

    if fastq_file[-3:] == '.gz':
        fastq_in = gzip.open(fastq_file, 'rt')
    elif fastq_file[-4:] == '.bz2':
        fastq_in = bz2.open(fastq_file, 'rt')
    else:
        fastq_in = open(fastq_file)

    header = fastq_in.readline()
    while header:
        seq = fastq_in.readline()
        mid = fastq_in.readline()
        qual = fastq_in.readline()

        if options.length_min is not None:
            if len(seq)-1 >= options.length_min:
                print('%s%s%s%s' % (header,seq,mid,qual), end='')

        header = fastq_in.readline()

    fastq_in.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
