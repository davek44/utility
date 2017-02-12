#!/usr/bin/env python
from optparse import OptionParser
import bz2
import gzip

'''
fastq_trim.py

Filter a FASTQ file for various properties, like read length.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <trim_length> <fastq_file>'
    parser = OptionParser(usage)
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

        # trim
        seq = seq[:trim_length]
        qual = qual[:trim_length]                  

        print('%s%s%s%s' % (header,seq,mid,qual), end='')

        header = fastq_in.readline()

    fastq_in.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
