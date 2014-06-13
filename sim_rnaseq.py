#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import lognorm, poisson
import pdb, random, subprocess, sys
import pysam
import gff, stats

################################################################################
# sim_rnaseq.py
#
# Simulate an RNA-Seq dataset from a GTF file and a FASTA file of transcript
# sequences.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf> <fasta>'
    parser = OptionParser(usage)
    parser.add_option('-b', dest='bam_length', help='Obtain read length via sampling a distribution from a BAM file [Default: %default]')
    parser.add_option('-e', dest='error_rate', type='float', default=0, help='Error rate (uniform on reads) [Default: %default]')
    parser.add_option('-f', dest='fpkm_file', help='Cufflinks .fpkm_tracking file to use for FPKMs [Default: %default]')
    parser.add_option('-l', dest='read_length', type='int', default=30, help='Read length [Default: %default]')
    parser.add_option('-n', dest='num_reads', type='int', default=100000, help='Number of reads [Default: %default]')
    parser.add_option('-o', dest='output_prefix', default='reads', help='Output files prefix [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide GTF file and fasta file')
    else:
        gtf_file = args[0]
        fasta_file = args[1]

    if options.bam_length:
        read_length_distribution = bam_length_distribution(options.bam_length)
    else:
        read_length_distribution = {options.read_length:1}

    # read GTF gene_id to transcript_id's mapping
    g2t = gff.g2t(gtf_file)

    # get transcript lengths
    transcript_lengths = {}
    for line in open(gtf_file):
        a = line.split('\t')
        if a[2] == 'exon':
            transcript_id = gff.gtf_kv(a[8])['transcript_id']
            transcript_lengths[transcript_id] = transcript_lengths.get(transcript_id,0) + int(a[4])-int(a[3])+1

    if options.fpkm_file:
        transcript_copies = {}
        fpkm_in = open(options.fpkm_file)
        line = fpkm_in.readline()
        for line in fpkm_in:
            a = line.split('\t')
            transcript_copies[a[0]] = float(a[9])
        fpkm_in.close()

        if sum(transcript_copies.values()) == 0:
            print >> sys.stderr, 'FPKM file shows no expression. Exiting.'
            exit(1)
    else:
        # sample gene copies
        gene_copies_raw = lognorm.rvs(1,size=len(g2t))
        gene_copies_raw_sum = sum(gene_copies_raw)
        gene_copies = dict(zip(g2t.keys(), [gcr/gene_copies_raw_sum for gcr in gene_copies_raw]))

        # sample transcript copies
        transcript_copies = {}
        for gene_id in g2t:
            relative_copies = dict(zip(g2t[gene_id], lognorm.rvs(1,size=len(g2t[gene_id]))))
            relative_sum = sum(relative_copies.values())
            for transcript_id in g2t[gene_id]:
                transcript_copies[transcript_id] = gene_copies[gene_id]*relative_copies[transcript_id]/relative_sum

    # determine transcript probabilities as a function of copy and length
    transcript_weights = {}
    for transcript_id in transcript_copies:
        if transcript_lengths[transcript_id] >= min(read_length_distribution.keys()):
            weight = 0
            for read_length in read_length_distribution:
                weight += read_length_distribution[read_length]*transcript_copies[transcript_id]*(transcript_lengths[transcript_id]-read_length+1)

            if weight > 0:
                transcript_weights[transcript_id] = weight
    weights_sum = sum(transcript_weights.values())
    transcript_probs = dict([(tid,transcript_weights[tid]/weights_sum) for tid in transcript_weights])

    # open fasta file
    fasta = pysam.Fastafile(fasta_file)

    # open output files
    fastq_out = open('%s.fastq' % options.output_prefix, 'w')
    gff_out = open('%s_txome.gff' % options.output_prefix, 'w')

    # for each transcript
    read_index = 1
    for transcript_id in transcript_probs:
        expected_reads = transcript_probs[transcript_id]*options.num_reads
        if expected_reads == 0:
            sampled_reads = 0
        else:
            sampled_reads = poisson.rvs(expected_reads)

        for s in range(sampled_reads):
            read_length = sample_read_length(read_length_distribution)
            if transcript_lengths[transcript_id] > read_length:
                pos = random.randint(0, transcript_lengths[transcript_id]-read_length)
                seq = fasta.fetch(transcript_id, pos, pos+read_length).upper()
                if seq:
                    eseq = inject_errors(seq, options.error_rate)

                    print >> fastq_out, '@read%d\n%s\n+\n%s' % (read_index,eseq,'I'*read_length)
                    print >> gff_out, '\t'.join([transcript_id, 'sim', 'read', str(pos+1), str(pos+read_length), '.', '+', '.', 'read%d'%read_index])

                    read_index += 1
                else:
                    print >> sys.stderr, 'Missing fasta sequence %s:%d-%d' % (transcript_id,pos,(pos+read_length))

    fastq_out.close()
    gff_out.close()

    # map back to genome
    subprocess.call('tgff_cgff.py -c %s %s_txome.gff > %s_genome.gff' % (gtf_file, options.output_prefix, options.output_prefix), shell=True)


################################################################################
# bam_length_distribution
#
# Input
#  bam_file: BAM file
#
# Output
#  length_probs:  Dict mapping read lengths to probabilities.
################################################################################
def bam_length_distribution(bam_file):
    length_counts = {}

    for aligned_read in pysam.Samfile(bam_file, 'rb'):
        length_counts[aligned_read.qlen] = length_counts.get(aligned_read.qlen,0) + 1.0/aligned_read.opt('NH')

    counts_sum = float(sum(length_counts.values()))
    length_probs = {}
    for rlen in length_counts:
        length_probs[rlen] = length_counts[rlen]/counts_sum

    return length_probs
    

################################################################################
# inject_errors
#
# 
################################################################################
def inject_errors(seq, error_rate):
    eseq = list(seq)
    for i in range(len(eseq)):
        if random.random() < error_rate:
            if eseq[i] == 'A':
                eseq[i] = random.choice(['C','G','T'])
            elif eseq[i] == 'C':
                eseq[i] = random.choice(['A','G','T'])
            elif eseq[i] == 'G':
                eseq[i] = random.choice(['A','C','T'])
            elif eseq[i] == 'T':
                eseq[i] = random.choice(['A','C','G'])
    return ''.join(eseq)


################################################################################
# sample_read_length
#
# Input
#  read_length_distribution:
#
# Output
#  read_length
################################################################################
def sample_read_length(read_length_distribution):
    read_lengths = read_length_distribution.keys()
    probs = [read_length_distribution[rl] for rl in read_lengths]
    read_length = stats.sample_probs(read_lengths, probs, 1)[0]
    return read_length


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
