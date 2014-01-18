#!/usr/bin/env python
from optparse import OptionParser
from glob import glob
import gzip, os, pdb
import gff, dna

################################################################################
# transcripts_fasta.py
#
# Make a fasta file of transcripts from the gtf file
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <genome fasta> <transcripts gtf>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide genome fasta file and transcripts gtf file')
    else:
        genome_fasta = args[0]
        transcripts_gtf = args[1]

    transcript_seqs = {}
    transcript_genes = {}

    if genome_fasta[-2:] == 'gz':
        genome_open = gzip.open(genome_fasta)
    else:
        genome_open = open(genome_fasta)

    # process chromosomes
    chrom = ''
    line = genome_open.readline()
    while line:
        if line[0] == '>':
            if chrom:
                process_chrom(transcripts_gtf, chrom, seq, transcript_seqs, transcript_genes)

            chrom = line[1:].rstrip()
            seq = ''
        else:
            seq += line.rstrip()
        line = genome_open.readline()
    process_chrom(transcripts_gtf, chrom, seq, transcript_seqs, transcript_genes)

    # print fasta
    for tid in transcript_seqs:
        print '>%s gene=%s\n%s' % (tid,transcript_genes[tid],transcript_seqs[tid])


################################################################################
# process_chrom
#
# Build up transcript_seqs and transcript_genes hashes for the chromosome
# specified.
################################################################################
def process_chrom(transcripts_gtf, chrom, seq, transcript_seqs, transcript_genes):
    # find chr transcripts
    for line in open(transcripts_gtf):
        a = line.split('\t')
        if a[0] == chrom:
            kv = gff.gtf_kv(a[8])
            tid = kv['transcript_id']
            gid = kv['gene_id']

            exon_start = int(a[3])
            exon_end = int(a[4])

            exon_seq = seq[exon_start-1:exon_end]
            if a[6] == '+':
                transcript_seqs[tid] = transcript_seqs.get(tid,'') + exon_seq
            else:
                transcript_seqs[tid] = dna.rc(exon_seq) + transcript_seqs.get(tid,'')

            transcript_genes[tid] = gid


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
