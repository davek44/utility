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
    parser.add_option('-a', dest='add_coords_header', default=False, action='store_true', help='Add BED coordinates to the FASTA header [Default: %default]')
    parser.add_option('-c', dest='chr_dir', default='', help='Directory of chromosome files named according to the first column of the gff file')
    parser.add_option('-f', dest='fasta_file', default='%s/research/common/data/genomes/hg19/sequence/hg19.fa' % os.environ['HOME'], help='Fasta file [Default: %default]')
    parser.add_option('-l', dest='length_match', default=None, type='int', help='Match all entries to be the same length [Default: %default]')
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
# Print sequence features for the given header and seq from the given bed file.
################################################################################
def header_bed(header, seq, bed_file, options):
    for line in open(bed_file):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        if a[0] == header:
            # determine start and end 
            feat_start = int(a[1])
            feat_end = int(a[2])
            if options.length_match:
                feat_mid = int(0.5*feat_start + 0.5*feat_end)
                feat_start = feat_mid - options.length_match/2
                feat_end = feat_mid + options.length_match/2

            # determine strand
            feat_strand = '+'
            if len(a) > 5 and  a[5] == '-':
                feat_strand = '-'

            # determine header
            if options.add_coords_header:
                feat_header = '%s:%s-%s:%s' % (header,a[1],a[2],feat_strand)
            else:
                feat_header = ''
                if len(a) > 3 and a[3] != '.':
                    feat_header = a[3]
            
            # determine sequence
            feat_seq = ''
            
            # if negative index, start with N's
            if feat_start < 0:
                feat_seq += 'N'*(-feat_start)
                feat_start = 0

            # grab the genome sequence
            feat_seq += seq[feat_start:feat_end]

            # if it's too short, extend with N's
            if options.length_match and len(feat_seq) < options.length_match:
                feat_seq += 'N'*(options.length_match - len(feat_seq))

            # reverse complement
            if feat_strand == '-':
                feat_seq = dna.rc(seq[feat_start:feat_end])

            #print '>%s\n%s' % (feat_header, feat_seq)
            print '>%s' % feat_header
            i = 0
            while i < len(feat_seq):
                print feat_seq[i:i+60]
                i += 60


################################################################################
# header_bed_id
#
# Print sequence features for the given header and seq from the given bed file,
# merging features with the same ID.
################################################################################
def header_bed_id(header, seq, bed_file, options):
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

    for head_id in header_seqs:
        print '>%s\n%s' % (head_id,header_seqs[head_id])


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
