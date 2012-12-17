#!/usr/bin/env python
from optparse import OptionParser
from pygr import worldbase
from pygr.sequence import *
import pdb, sys
import dna, gff

################################################################################
# transmap_lncrna_mm9.py
#
# TransMap a specified lncRNA transcript to all genomes.
# That is,
#  1. Grab the lncRNA exon coordinates.
#  2. Search them against the 30 way multiz alignment via Pygr
#  3. For each genome, sort and merge those alignments
#  4. For each genome, if there is enough, spit out a gff and fasta file
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <transcript id>'
    parser = OptionParser(usage)
    parser.add_option('-a', dest='align_t', type='float', default=0.01, help='Minimum % of the transcript that must align [Default: %default]')
    parser.add_option('-l', dest='lncrna_gtf', default='/Users/dk/research/common/data/lncrna_mm9/lnc_catalog.gtf', help='lncRNA gtf file [Default: %default]')
    parser.add_option('-m', dest='merge_t', type='int', default=40, help='Minimum distance between alignment blocks to merge into a single exon [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide transcript id')
    else:
        transcript_id = args[0]

    # get human genome
    mm9 = worldbase.Bio.Seq.Genome.MOUSE.mm9()

    # get gene exon intervals
    gene_ivals = []
    transcript_length = 0
    for line in open(options.lncrna_gtf):
        a = line.split('\t')
        kv = gff.gtf_kv(a[8])
        if kv['transcript_id'] == transcript_id:
            chrom = a[0]
            start = int(a[3])
            end = int(a[4])
            strand = 1*(a[6]=='+') - 1*(a[6]=='-') # assuming all orientations are the same
            gene_id = kv['gene_id']

            gene_ivals.append(mm9[chrom][start:end])
            transcript_length += end-start

    # get mm9 msa
    msa = worldbase.Bio.MSA.UCSC.mm9_multiz30way()

    # map returned sequences back to genome name
    idDict = ~(msa.seqDict)

    # hash alignments by genome
    genome_blocks = {}
    for gi in gene_ivals:
        for src, dest, edg in msa[gi].edges():
            genome_blocks.setdefault(idDict[dest],[]).append(dest)
            #print repr(gi), repr(src), repr(dest), idDict[dest], edg.length()

    # check for enough alignment
    for gen_chr in genome_blocks.keys():
        aligned_nt = sum([b.stop-b.start for b in genome_blocks[gen_chr]])
        print gen_chr, aligned_nt, float(aligned_nt)/transcript_length
        if aligned_nt < options.align_t*transcript_length:
            del genome_blocks[gen_chr]

    # for each genome
    worldbase_genomes = worldbase.dir('Bio.Seq.Genome')
    for gen_chr in genome_blocks:
        genome_blocks[gen_chr].sort(block_cmp)

        # make gtf lines / merge alignments
        b = genome_blocks[gen_chr][0]
        gff_strand = '+'*(b.orientation==strand) + '-'*(b.orientation!=strand)
        gff_cols = [[b.id, 'PygrTransMap', 'exon', b._abs_interval[0]+1, b._abs_interval[1], '.', gff_strand, '.', 'gene_id "%s"; transcript_id "%s"; exon_number "1";' % (gene_id, transcript_id,)]]
        exon_num = 2
        for i in range(1,len(genome_blocks[gen_chr])):
            if gff_cols[-1][4] + options.merge_t >= genome_blocks[gen_chr][i]._abs_interval[0]:
                # merge with prior
                gff_cols[-1][4] = genome_blocks[gen_chr][i]._abs_interval[1]
            else:
                # add new exon
                b = genome_blocks[gen_chr][i]
                gff_cols.append([b.id, 'PygrTransMap', 'exon', b._abs_interval[0]+1, b._abs_interval[1], '.', gff_strand, '.', 'gene_id "%s"; transcript_id "%s"; exon_number "%d";' % (gene_id, transcript_id, exon_num)])
                exon_num += 1

        # print gtf
        gtf_out = open('%s_%s.gtf' % (transcript_id, gen_chr), 'w')
        for gc in gff_cols:
            print >> gtf_out, '\t'.join([str(c) for c in gc])
        gtf_out.close()

        # get genomic sequence
        gen = gen_chr[:gen_chr.find('.')]
        chrom = gen_chr[gen_chr.find('.')+1:]
        wb_gen = [wgen for wgen in worldbase_genomes if wgen.find(gen) != -1]
        if len(wb_gen) > 1:
            print >> sys.stderr, 'Detected >1 worldbase genome matching %s' % gen
            print >> sys.stderr, ' '.join(wb_gen)         
        gen_seq = worldbase.__call__(wb_gen[0])

        # get transcript sequence
        seq = ''
        for gc in gff_cols:
            seq += str(gen_seq[gc[0]][gc[3]-1:gc[4]])
        if gff_cols[0][6] == '-':
            seq = dna.rc(seq)

        # print fasta
        fasta_out = open('%s_%s.fa' % (transcript_id, gen_chr), 'w')
        print >> fasta_out, '>%s_gene=%s_%s\n%s' % (transcript_id,gene_id,gen_chr,seq)
        fasta_out.close()


def block_cmp(x, y):
    return x._abs_interval[0] - y._abs_interval[0]


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
