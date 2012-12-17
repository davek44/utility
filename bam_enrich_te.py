#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import binom
import gzip, os, random, subprocess, sys
import pysam
import gff

################################################################################
# te_bam_enrich.py
#
# Compute the enrichment of aligned reads in a BAM file in transposable
# element families.
################################################################################

home_dir = '/Users/dk' # /n/home03/dkelley

linc_gtf = '%s/research/common/data/lncrna/linc_catalog_trans.gtf' % home_dir

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam file>'
    parser = OptionParser(usage)
    parser.add_option('-f', dest='feature', help='Use only rows from the feature gff where the 3rd column matches feature')
    parser.add_option('-g', dest='genome_fa', default='%s/research/common/data/genomes/hg19/sequence/hg19.fa' % home_dir, help='Genome fasta file [Default: %default]')
    parser.add_option('-l', dest='lincrna', default=False, action='store_true', help='Only use lincRNA promoter-associated TEs [Default: %default]')
    parser.add_option('-n', dest='no_nh', default=False, action='store_true', help='BAM alignments lack the NH tag for multi-mapping reads [Default: %default]')
    parser.add_option('-r', dest='repeats_gff', default='%s/research/common/data/genomes/hg19/annotation/repeatmasker/hg19.fa.out.tp.gff' % home_dir)
    (options,args) = parser.parse_args()

    if len(args) < 1:
        parser.error('Must provide a gff file for the feature of interest.')
    else:
        bam_files = args

    # count genomic bp
    genome_bp = count_fasta(options.genome_fa)

    if options.lincrna:
        linc_prom_gtf = '%s_prom.gtf' % os.path.splitext(linc_gtf)[0]
        gff.promoters(linc_gtf, 2000, 200, output_file=linc_prom_gtf)

        tmp_gff = 'te_linc_%d.gff' % random.randint(0,1000000)
        p = subprocess.Popen('intersectBed -wa -u -a %s -b %s > %s' % (options.repeats_gff, linc_prom_gtf, tmp_gff), shell=True)
        os.waitpid(p.pid,0)
        options.repeats_gff = tmp_gff

    # hash counted repeat genomic bp
    te_lengths = measure_te(options.repeats_gff)

    te_counts = {}
    num_aligned_reads = 0
    for bam_file in bam_files:
        # count # aligned reads
        # hash multi-mapping reads
        multi_reads = {}
        if options.no_nh:
            bam_in = pysam.Samfile(bam_file, 'rb')
            for read in bam_in:
                multi_reads[read.qname] = multi_reads.get(read.qname,0) + 1
            bam_in.close()

            num_aligned_reads += len(multi_reads)
            for qname in multi_reads.keys():
                if multi_reads[qname] == 1:
                    del multi_reads[qname]

        else:
            num_aligns = 0
            bam_in = pysam.Samfile(bam_file, 'rb')
            for read in bam_in:
                num_aligns += 1
                if read.opt('NH') > 1:
                    multi_reads[read.qname] = read.opt('NH')
            bam_in.close()
            num_aligned_reads += num_aligns - sum(multi_reads.values()) + len(multi_reads)

        # intersect (require 50% of read)
        proc = subprocess.Popen('intersectBed -wo -bed -f 0.5 -abam %s -b %s' % (bam_file,options.repeats_gff), shell=True, stdout=subprocess.PIPE)

        # hash read counts by TE family
        line = proc.stdout.readline()
        while line:
            a = line.split('\t')
            te_kv = gff.gtf_kv(a[14])

            if not a[3] in multi_reads:
                read_inc = 1.0
            else:
                read_inc = 1.0/multi_reads[a[3]]

            te_counts[(te_kv['repeat'],te_kv['family'])] = te_counts.get((te_kv['repeat'],te_kv['family']),0.0) + read_inc
            te_counts[('*',te_kv['family'])] = te_counts.get(('*',te_kv['family']),0.0) + read_inc
            te_counts[('*','*')] = te_counts.get(('*','*'),0.0) + read_inc

            line = proc.stdout.readline()
        proc.communicate()

    if options.lincrna:
        os.remove(tmp_gff)

    # compute stats, print table
    for (rep,fam) in te_counts:
        te_p = float(te_lengths[(rep,fam)]) / genome_bp

        if te_counts[(rep,fam)] > te_p*num_aligned_reads:
            p_val = binom.sf(int(te_counts[(rep,fam)])-1, num_aligned_reads, te_p)
        else:
            p_val = binom.cdf(int(te_counts[(rep,fam)]), num_aligned_reads, te_p)

        if te_p*num_aligned_reads > 0:
            fold = te_counts[(rep,fam)]/(te_p*num_aligned_reads)
        else:
            fold = 0

        cols = (rep, fam, te_lengths[(rep,fam)], te_counts[(rep,fam)], te_p, fold, p_val)

        print '%-18s %-18s %10d %10d %10.2e % 10.3f %10.2e' % cols


################################################################################
# count_fasta
#
# Count the number of non-N nucleotides in a fasta file
################################################################################
def count_fasta(fasta_file):
    if fasta_file[-2:] == 'gz':
        f = gzip.open(fasta_file)
    else:
        f = open(fasta_file)

    line = f.readline()
    bp = 0
    while line:
        if line[0] != '>':
            bp += len([nt for nt in line.rstrip() if nt != 'N'])
        line = f.readline()
    return bp


################################################################################
# measure_te
#
# Hash the number of bp covered by various repeats in the RepeatMasker gff file
# and the lincRNA gtf file.
################################################################################
def measure_te(rm_file):
    repeat_bp = {}
    for line in open(rm_file):
        a = line.split('\t')

        kv = gff.gtf_kv(a[8])
        rep = kv['repeat']
        family = kv['family']

        length = int(a[4]) - int(a[3]) + 1

        repeat_bp[(rep,family)] = repeat_bp.get((rep,family),0) + length
        repeat_bp[('*',family)] = repeat_bp.get(('*',family),0) + length
        repeat_bp[('*','*')] = repeat_bp.get(('*','*'),0) + length

    return repeat_bp


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
