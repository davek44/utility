#!/usr/bin/env python
from optparse import OptionParser
import gzip, os, pdb, re, sys
import dna, gff

################################################################################
# gtf2prom.py
#
# Produce a GFF file and a fasta file corresponding to the promoter of the
# genes in a gtf fle.
################################################################################

hg19_fa = '/Users/dk/research/common/data/genomes/hg19/sequence/hg19.fa'

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf file>'
    parser = OptionParser(usage)
    parser.add_option('-l', dest='length', type='int', default=2000, help='Promoter length [Default: %default]')
    parser.add_option('-n', dest='acgt_t', type='float', default=0.9, help='Proportion of non-N nt\'s allowed [Default: %default]')
    parser.add_option('-o', dest='output_pre', help='Output file prefix')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide gtf file')
    else:
        gtf_file = args[0]

    promoters = get_promoters(gtf_file, options.length)
    output_promoters(promoters, options.length, options.acgt_t, options.output_pre)


################################################################################
# acgt_pct
################################################################################
def acgt_pct(seq):
    acgt = 0
    for nt in seq:
        if 'ACGTacgt'.find(nt) != -1:
            acgt += 1
    return float(acgt)/len(seq)


################################################################################
# find_promoter
#
# Determine the position of the promoter from the first exon of the gene
################################################################################
def find_promoter(gene_id, exons, promoter_length):
    strand = exons[0].split()[6]
    exon_starts = [x.split()[3] for x in exons]

    if strand == '+':
        exon_starts.sort()
    else:
        exon_starts.sort(reverse=True) 

    initial_exon = [x for x in exons if x.split()[3] == exon_starts[0]][0]

    return Promoter(initial_exon, promoter_length)


################################################################################
# get_promoters
#
# Make a list of promoters from the linc catalog
################################################################################
def get_promoters(gtf_file, promoter_length):
    promoters = []

    gene_id = ''
    for line in open(gtf_file):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        this_gene_id = gff.gtf_kv(a[8])['gene_id']
        if this_gene_id != gene_id:
            if gene_id:
                promoters.append(find_promoter(gene_id, exons, promoter_length))
            gene_id = this_gene_id                
            exons = [line]
        else:
            exons.append(line)

    promoters.append(find_promoter(gene_id, exons, promoter_length))

    return promoters


################################################################################
# output_promoters
#
# Output fasta and gff files for the promoters.
################################################################################
def output_promoters(promoters, promoter_length, acgt_t, output_pre):
    # hash promoters by chr
    chr_hash = {}
    for prom in promoters:
        chr_hash.setdefault(prom.chr, []).append(prom)

    if output_pre:
        out_fa = open('%s.fa' % output_pre,'w')
        out_gff = open('%s.gff' % output_pre,'w')
    else:
        out_fa = open('promoters.fa','w')
        out_gff = open('promoters.gff','w')

    if os.path.isfile(hg19_fa):        
        genome_in = open(hg19_fa)
    elif os.path.isfile(hg19_fa+'.gz'):
        genome_in = gzip.open(hg19_fa+'.gz')
    else:
        print >> sys.stderr, 'No genome %s' % hg19_fa
        exit(1)

    chrom = ''
    line = genome_in.readline()
    while line:
        if line[0] == '>':
            if chrom:
                process_chr(chrom, chr_seq, chr_hash.get(chrom,[]), out_fa, out_gff, promoter_length, acgt_t)

            chrom = line[1:].rstrip()
            chr_seq = ''
        else:
            chr_seq += line.rstrip()
        line = genome_in.readline()
    process_chr(chrom, chr_seq, chr_hash.get(chrom,[]), out_fa, out_gff, promoter_length, acgt_t)

    out_fa.close()
    out_gff.close()


################################################################################
# process_chr
#
# Print promoter information for a single chromosome.
################################################################################
def process_chr(chrom, seq, promoters, out_fa, out_gff, promoter_length, acgt_t):
    # grab promoters
    for prom in promoters:
        if prom.strand == '+':
            prom_seq = seq[prom.start:prom.start+promoter_length]
        else:
            prom_seq = dna.rc(seq[prom.start:prom.start+promoter_length])
        if acgt_pct(prom_seq) > acgt_t:
            print >> out_fa, '>%s\n%s' % (prom.gtf_kv['transcript_id'], prom_seq)
            gff_dat = [chrom, '.', 'promoter', str(prom.start+1), str(prom.start+promoter_length+1-1), '.', prom.strand, '.', gff.kv_gtf(prom.gtf_kv)]
            print >> out_gff, '\t'.join(gff_dat)


################################################################################
# Promoter
################################################################################
class Promoter:
    def __init__(self, exon_gtf, promoter_length):
        a = exon_gtf.split('\t')
        a[-1] = a[-1].rstrip()

        self.gtf_kv = gff.gtf_kv(a[8])
        self.chr = a[0]
        self.strand = a[6]

        if self.strand == '+':
            self.start = max(0, int(a[3]) - promoter_length)
        else:
            self.start = int(a[4])


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
