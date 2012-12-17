#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import binom
import gzip, os, random, subprocess, sys, tempfile
import pysam
import gff

################################################################################
# te_bam_enrich.py
#
# Compute the enrichment of aligned reads in a BAM file in transposable
# element families relative to a null model determined by a control BAM file.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam file> <control bam file>'
    parser = OptionParser(usage)

    parser.add_option('-g', dest='gff_file', default=None, help='Filter the TEs by overlap with genes in the given gff file [Default: %default]')

    parser.add_option('-n', dest='no_nh', default=False, action='store_true', help='BAM alignments lack the NH tag for multi-mapping reads [Default: %default]')

    parser.add_option('-r', dest='repeats_gff', default='%s/research/common/data/genomes/hg19/annotation/repeatmasker/hg19.fa.out.tp.gff' % os.environ['HOME'])
    (options,args) = parser.parse_args()

    if len(args) < 2:
        parser.error(usage)
    else:
        bam_file = args[0]
        control_bam_file = args[1]

    # filter TEs by gff file
    if options.gff_file:
        te_gff_fd, te_gff_file = tempfile.mkstemp()

        p = subprocess.Popen('intersectBed -u -a %s -b %s > %s' % (options.repeats_gff, options.gff_file, te_gff_file), shell=True)
        os.waitpid(p.pid,0)
        options.repeats_gff = te_gff_file

    # count read alignments to TEs
    te_counts, num_aligned_reads = count_bam(bam_file, options.repeats_gff, options.no_nh)
    te_counts_control, num_aligned_reads_control = count_bam(control_bam_file, options.repeats_gff, options.no_nh)

    # compute stats, print table
    for (rep,fam) in te_counts:
        if not te_counts_control.has_key((rep,fam)):
            cols = (rep, fam, te_counts[(rep,fam)], '-', '-', '-')
            print '%-18s %-18s %10d %10s %10s %10s' % cols
        else:
            te_p = float(te_counts_control[(rep,fam)]) / num_aligned_reads_control

            if te_counts[(rep,fam)] > te_p*num_aligned_reads:
                p_val = binom.sf(int(te_counts[(rep,fam)])-1, num_aligned_reads, te_p)
            else:
                p_val = binom.cdf(int(te_counts[(rep,fam)]), num_aligned_reads, te_p)

            if te_p*num_aligned_reads > 0:
                fold = te_counts[(rep,fam)]/(te_p*num_aligned_reads)
            else:
                fold = 0

            cols = (rep, fam, te_counts[(rep,fam)], te_p, fold, p_val)
            print '%-18s %-18s %10d %10.2e %10.3f %10.2e' % cols


################################################################################
# count_bam
#
# Hash read alignment counts by TEs, and report the overall number of aligned
# reads.
################################################################################
def count_bam(bam_file, repeats_gff, no_nh):
    te_counts = {}
    num_aligned_reads = 0

    # count # aligned reads
    # hash multi-mapping reads
    multi_reads = {}
    if no_nh:
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
    proc = subprocess.Popen('intersectBed -wo -bed -f 0.5 -abam %s -b %s' % (bam_file,repeats_gff), shell=True, stdout=subprocess.PIPE)

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

    return te_counts, num_aligned_reads


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
