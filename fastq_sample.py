#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
import gzip
import random

################################################################################
# fastq_sample.py
#
#
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <sample_num> <in1_fq> <in2_fq> <out1_fq> <out2_fq>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 5:
        parser.error('Must provide sample number, two input FASTQs, and two output FASTQs')
    else:
        sample_num = int(args[0])
        in1_fq = args[1]
        in2_fq = args[2]
        out1_fq = args[3]
        out2_fq = args[4]

    ##################################################
    # count fragments

    fragment_num = 0
    in1_open = open_maygz(in1_fq)
    header = in1_open.readline()
    while header:
        seq = in1_open.readline()
        mid = in1_open.readline()
        qual = in1_open.readline()
        fragment_num += 1
        header = in1_open.readline()
    in1_open.close()


    ##################################################
    # sample the indexes of reads to keep

    if sample_num < fragment_num:
        sampled_indexes = sorted(random.sample(range(fragment_num), sample_num))
    else:
        sampled_indexes = range(fragment_num)

    ##################################################
    # pass through files and output those reads

    in1_open = open_maygz(in1_fq)
    in2_open = open_maygz(in2_fq)
    out1_open = open(out1_fq, 'w')
    out2_open = open(out2_fq, 'w')

    fi = 0
    si = 0

    # read first entries
    header1, seq1, mid1, qual1 = read_fastq(in1_open)
    header2, seq2, mid2, qual2 = read_fastq(in2_open)

    while header1:
        # if selected
        if si < len(sampled_indexes) and fi == sampled_indexes[si]:
            # print reads
            print_read(out1_open, header1, seq1, mid1, qual1)
            print_read(out2_open, header2, seq2, mid2, qual2)

            # advance sampled index
            si += 1

        # advance fastq index
        fi += 1

        # read next entries
        header1, seq1, mid1, qual1 = read_fastq(in1_open)
        header2, seq2, mid2, qual2 = read_fastq(in2_open)

    in1_open.close()
    in2_open.close()
    out1_open.close()
    out2_open.close()


def open_maygz(input_file):
    ''' Open the file, which may be gzipped. '''
    if input_file[-3:] == '.gz':
        input_open = gzip.open(input_file)
    else:
        input_open = open(input_file)
    return input_open


def print_read(out_open, header, seq, mid, qual):
    print(header, file=out_open, end='')
    print(seq, file=out_open, end='')
    print(mid, file=out_open, end='')
    print(qual, file=out_open, end='')


def read_fastq(input_open):
    ''' Read a FASTQ entry from an open, possible gzipped file. '''

    header = input_open.readline()
    if header:
        seq = input_open.readline()
        mid = input_open.readline()
        qual = input_open.readline()

        if type(header) == bytes:
            header = header.decode('UTF-8')
            seq = seq.decode('UTF-8')
            mid = mid.decode('UTF-8')
            qual = qual.decode('UTF-8')
    else:
        seq, mid, qual = None, None, None

    return header, seq, mid, qual


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
