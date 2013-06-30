#!/usr/bin/env python
from optparse import OptionParser
import gzip, os, sys
import dna

################################################################################
# bed2fa.py
#
# Given a bed file and fasta file or chromosome fasta file directory, produce
# a fasta file of the bed entries.
#
# WARNING: Simple first attempt- no block support
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bed file>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='chr_dir', default='', help='Directory of chromosome files named according to the first column of the gff file')
    parser.add_option('-f', dest='fasta_file', default='/n/rinn_data1/indexes/human/hg19/bowtie2/hg19.fa', help='Fasta file [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error(usage)
    else:
        bed_file = args[0]

    if options.fasta_file:
        fasta_files = [options.fasta_file]
    elif options.chr_dir:
        bed_chrs = set([line.split()[0] for line in open(bed_file)])
        fasta_files = []
        for chrom in bed_chrs:
            fasta_files += glob.glob(options.chr_dir+'/%s*' % chrom)
    else:
        parser.error('Must provide fasta source')

    for fasta_file in fasta_files:
        if fasta_file[-3:] == '.gz':
            fasta_open = gzip.open(fasta_file)
        else:
            fasta_open = open(fasta_file)

        header = ''
        line = fasta_open.readline()
        while line:
            if line[0] == '>':
                if header:
                    header_bed(header, seq, bed_file, options)
                header = line[1:].split()[0]
                seq = ''
            else:
                seq += line.rstrip()
            line = fasta_open.readline()
        header_bed(header, seq, bed_file, options)


################################################################################
# header_bed
#
# Print sequence features for the given header and seq from the given bed file
################################################################################
def header_bed(header, seq, bed_file, options):
    header_seqs = {}
    for line in open(bed_file):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        if a[0] == header:
            feat_start = int(a[1])
            feat_end = int(a[2])
            head_id = a[3]
            
            feat_seq = seq[feat_start:feat_end]

            if a[5] == '+':
                header_seqs[head_id] = header_seqs.get(head_id,'') + feat_seq
            else:
                header_seqs[head_id] = dna.rc(feat_seq) + header_seqs.get(head_id,'')

    for header in header_seqs:
        print '>%s\n%s' % (header,header_seqs[header])


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
