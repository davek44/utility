#!/usr/bin/env python
from optparse import OptionParser
import gzip, sys

############################################################
# sum_fasta
#
# Sum the lengths of all sequences listed in a multi fasta
# file.
############################################################

############################################################
# main
############################################################
def main():
    parser = OptionParser()
    parser.add_option('-i', action='store_true', dest='sum_indiv', default=False, help='Sum all entries in fasta file individually [Default: %default]')
    parser.add_option('-n', action='store_true', dest='ignore_n', default=False, help='Ignore N basepairs [Default: %default]')
    parser.add_option('-q', action='store_true', dest='fastq', default=False, help='File is fastq [Default: %default]')
    parser.add_option('-r', action='store_true', dest='raw_only', default=False, help='Print only the raw count rather than a string version, too [Default: %default]')
    (options,args) = parser.parse_args()

    if options.sum_indiv:
        for faf in args:
            sums_fasta(faf, options.ignore_n, options.raw_only)

    else:
        for faf in args:
            print('%s\t' % faf, end='')
            if options.fastq:
                sum_fastq(faf, True, options.ignore_n, options.raw_only)
            else:
                sum_fasta(faf, True, options.ignore_n, options.raw_only)


############################################################
# sum_fasta
############################################################
def sum_fasta(fasta_file, do_print=False, ignore_n=False, raw_only=False):
    if fasta_file[-3:] == '.gz':
        fasta_in = gzip.open(fasta_file)
    else:
        fasta_in = open(fasta_file)

    fsum = 0
    line = fasta_in.readline()
    while line:
        if line[0] != '>':
            if ignore_n:
                fsum += len([nt for nt in line.rstrip() if nt != 'N'])
            else:
                fsum += len(line.rstrip())
        line = fasta_in.readline()

    if do_print:
        output_sum(fsum, raw_only)

    return fsum


############################################################
# sum_fastq
############################################################
def sum_fastq(fastq_file, do_print=False, ignore_n=False, raw_only=False):
    if fastq_file[-3:] == '.gz':
        fastq_in = gzip.open(fastq_file)
    else:
        fastq_in = open(fastq_file)

    # sum sequences
    fsum = 0
    header = fastq_in.readline()
    while header:
        seq = fastq_in.readline().rstrip()
        garbage = fastq_in.readline()
        qual = fastq_in.readline()

        if ignore_n:
            fsum += len([nt for nt in seq if nt != 'N'])
        else:
            fsum += len(seq)

        header = fastq_in.readline()

    if do_print:
        output_sum(fsum, raw_only)

    return fsum

############################################################
# sums_fasta
############################################################
def sums_fasta(fasta_file, ignore_n=False, raw_only=False):
    if fasta_file[-3:] == '.gz':
        fasta_in = gzip.open(fasta_file)
    else:
        fasta_in = open(fasta_file)

    fsum = 0
    line = fasta_in.readline()
    while line:
        if line[0] == '>':
            if fsum:
                output_sum(fsum, raw_only)
            fsum = 0
            print(line.rstrip()+'\t', end='')
        else:
            if ignore_n:
                fsum += len([nt for nt in line.rstrip() if nt != 'N'])
            else:
                fsum += len(line.rstrip())
        line = fasta_in.readline()

    if fsum:
        output_sum(fsum, raw_only)


############################################################
# output_sum
############################################################
def output_sum(sum, raw_only):
    # output
    print(sum, end='')
    if raw_only:
        print('')
    else:
        if sum > 1000000000:
            print('(%.4f Gb)' % (sum/1000000000.0))
        elif sum > 1000000:
            print('(%.4f Mb)' % (sum/1000000.0))
        elif sum > 1000:
            print('(%.4f Kb)' % (sum/1000.0))
        else:
            print('')


############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
