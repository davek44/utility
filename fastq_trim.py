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

    if len(args) != 2:
        parser.error('Must provide trim length and FASTQ file')
    else:
        trim_length = int(args[0])
        fastq_file = args[1]

    if fastq_file[-3:] == '.gz':
        fastq_in = gzip.open(fastq_file, 'rt')
    elif fastq_file[-4:] == '.bz2':
        fastq_in = bz2.open(fastq_file, 'rt')
    else:
        fastq_in = open(fastq_file)

    header = fastq_in.readline().rstrip()
    while header:
        seq = fastq_in.readline().rstrip()
        mid = fastq_in.readline().rstrip()
        qual = fastq_in.readline().rstrip()

        # trim
        seq = seq[:trim_length]
        qual = qual[:trim_length]                  

        print('%s\n%s\n%s\n%s' % (header,seq,mid,qual))

        header = fastq_in.readline().rstrip()

    fastq_in.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
