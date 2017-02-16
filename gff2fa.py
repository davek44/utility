#!/usr/bin/env python
from optparse import OptionParser
import gzip, os, pdb, sys
import dna, gff

################################################################################
# gff2fa.py
#
# Given a gff file and fasta file or chromosome fasta file directory, produce
# a fasta file of the gff entries.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gff file>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='chr_dir', default='', help='Directory of chromosome files named according to the first column of the gff file')
    parser.add_option('-f', dest='fasta_file', default='%s/research/common/data/genomes/hg19/sequence/hg19.fa' % os.environ['HOME'], help='Fasta file with entries named according to the first column of the gff file [Default: %default]')
    parser.add_option('-g', dest='gene_too', default=False, action='store_true', help='Print transcript id and gene id [Default: %default]')
    #parser.add_option('--gtf', dest='gtf', action='store_true', default=False, help='Input file is gtf and linked entries should be combined into a single sequence [Default: %default]')
    parser.add_option('--head', dest='header_key', default=None, help='GFF key to be used to merge GFF entries and label the fasta headers [Default: %default]')
    parser.add_option('-s', dest='split_lines', default=None, action='store_true', help='Split sequence across multiple lines [Default: %default]')
    parser.add_option('-x', dest='exon', action='store_true', default=False, help='Only include exon rows [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide gff file')
    else:
        gff_file = args[0]

    if options.fasta_file:
        fasta_files = [options.fasta_file]
    elif options.chr_dir:
        gff_chrs = set([line.split()[0] for line in open(gff_file)])
        fasta_files = []
        for chrom in gff_chrs:
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
                    header_gff(header, seq, gff_file, options)
                header = line[1:].split()[0]
                seq = ''
            else:
                seq += line.rstrip()
            line = fasta_open.readline()
        header_gff(header, seq, gff_file, options)


################################################################################
# header_gff
#
# Print sequence features for the given header and seq from the given gff file
################################################################################
def header_gff(header, seq, gff_file, options):
    header_seqs = {}
    for line in open(gff_file):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()
        if (not options.exon or a[2] == 'exon') and a[0] == header:
            try:
                kv = gff.gtf_kv(a[8])
            except:
                kv = {}

            head_id = kv.get(options.header_key,a[0]+':'+a[3]+'-'+a[4])
            #head_id = kv.get(options.header_key,a[8])

            if options.gene_too:
                head_id += ' gene=%s' % kv.get('gene_id','')

            feat_start = int(a[3])
            feat_end = int(a[4])

            feat_seq = seq[feat_start-1:feat_end]
            if a[6] == '+':
                header_seqs[head_id] = header_seqs.get(head_id,'') + feat_seq
            else:
                header_seqs[head_id] = dna.rc(feat_seq) + header_seqs.get(head_id,'')

    for header in header_seqs:
        print '>%s' % header
        if options.split_lines:
            i = 0
            while i < len(header_seqs[header]):
                print header_seqs[header][i:i+60]
                i += 60
        else:
            print header_seqs[header]


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
