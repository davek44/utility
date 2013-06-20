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
# element families.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam file>'
    parser = OptionParser(usage)

    parser.add_option('-g', dest='gff_file', default=None, help='Filter the TEs by overlap with genes in the given gff file [Default: %default]')

    parser.add_option('-r', dest='repeats_gff', default='%s/research/common/data/genomes/hg19/annotation/repeatmasker/hg19.fa.out.tp.gff' % os.environ['HOME'])
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide a gff file for the feature of interest.')
    else:
        bam_file = args[0]

    # compute size of search space
    if options.gff_file:
        genome_bp = count_gff(options.gff_file)
    else:
        genome_bp = count_hg19()

    # filter TEs and read alignments by gff file
    if options.gff_file:
        te_gff_fd, te_gff_file = tempfile.mkstemp(dir='%s/research/scratch' % os.environ['HOME'])

        subprocess.call('intersectBed -u -a %s -b %s > %s' % (options.repeats_gff, options.gff_file, te_gff_file), shell=True)
        options.repeats_gff = te_gff_file

        bam_gff_fd, bam_gff_file = tempfile.mkstemp(dir='%s/research/scratch' % os.environ['HOME'])
        subprocess.call('intersectBed -abam %s -b %s > %s' % (bam_file, options.gff_file, bam_gff_file), shell=True)
        bam_file = bam_gff_file

    # filter BAM file for mapping quality
    bam_mapq_fd, bam_mapq_file = tempfile.mkstemp(dir='%s/research/scratch' % os.environ['HOME'])
    bam_in = pysam.Samfile(bam_file, 'rb')
    bam_mapq_out = pysam.Samfile(bam_mapq_file, 'wb', template=bam_in)
    for aligned_read in bam_in:
        if aligned_read.mapq > 0:
            bam_mapq_out.write(aligned_read)
    bam_mapq_out.close()

    # hash counted repeat genomic bp
    te_lengths = measure_te(options.repeats_gff)

    # count fragments and hash multi-mappers
    num_fragments = 0
    multi_maps = {}
    paired_poll = {False:0, True:0}
    for aligned_read in pysam.Samfile(bam_mapq_file, 'rb'):
        if aligned_read.is_paired:
            num_fragments += 0.5/aligned_read.opt('NH')
        else:
            num_fragments += 1.0/aligned_read.opt('NH')

        if aligned_read.opt('NH') > 0:
            multi_maps[aligned_read.qname] = aligned_read.opt('NH')

        paired_poll[aligned_read.is_paired] += 1

    # guess paired-ness
    if paired_poll[True] > 0 and paired_poll[False] > 0:
        print >> sys.stderr, 'Paired-ness of the reads is ambiguous'
    if paired_poll[True] > paired_poll[False]:
        is_paired = True
    else:
        is_paired = False

    # hash read counts by TE family
    te_counts = {}
    proc = subprocess.Popen('intersectBed -wo -bed -abam %s -b %s' % (bam_mapq_file,options.repeats_gff), shell=True, stdout=subprocess.PIPE)
    for line in proc.stdout:
        a = line.split('\t')
        te_kv = gff.gtf_kv(a[14])

        if is_paired:
            read_inc = 0.5/multi_maps.get(a[3],1.0)
        else:
            read_inc = 1.0/multi_maps.get(a[3],1.0)

        te_counts[(te_kv['repeat'],te_kv['family'])] = te_counts.get((te_kv['repeat'],te_kv['family']),0.0) + read_inc
        te_counts[('*',te_kv['family'])] = te_counts.get(('*',te_kv['family']),0.0) + read_inc
        te_counts[('*','*')] = te_counts.get(('*','*'),0.0) + read_inc

    proc.communicate()

    # compute stats, print table
    for (rep,fam) in te_counts:
        te_p = float(te_lengths[(rep,fam)]) / genome_bp

        if te_counts[(rep,fam)] > te_p*num_fragments:
            p_val = binom.sf(int(te_counts[(rep,fam)])-1, int(num_fragments), te_p)
        else:
            p_val = binom.cdf(int(te_counts[(rep,fam)]), int(num_fragments), te_p)

        if te_p*num_fragments > 0:
            fold = te_counts[(rep,fam)]/(te_p*num_fragments)
        else:
            fold = 0

        cols = (rep, fam, te_lengths[(rep,fam)], te_counts[(rep,fam)], te_p, fold, p_val)

        print '%-18s %-18s %10d %10d %10.2e %10.3f %10.2e' % cols

    # clean
    os.remove(bam_mapq_file)
    if options.gff_file:
        os.remove(te_gff_file)
        os.remove(bam_gff_file)


################################################################################
# count_gff
#
# Count the number of bp in the limiting GFF file.
################################################################################
def count_gff(gff_file):
    gff_bp = 0
    p = subprocess.Popen('mergeBed -i %s' % gff_file, shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        a = line.split()
        gff_bp += int(a[2]) - int(a[1])
    p.communicate()

    return gff_bp


################################################################################
# count_hg19
#
# Count the number of bp in hg19 where TEs could be.
################################################################################
def count_hg19():
    chrom_sizes_file = '%s/research/common/data/genomes/hg19/assembly/human.hg19.genome' % os.environ['HOME']
    gap_bed_file = '%s/research/common/data/genomes/hg19/assembly/hg19_gaps.bed' % os.environ['HOME']
    valid_chrs = ['chr%d' % c for c in range(1,23)] + ['chrX','chrY']

    genome_bp = 0
    for line in open(chrom_sizes_file):        
        a = line.split()
        if len(a) > 0 and a[0] in valid_chrs:
            genome_bp += int(a[1])

    for line in open(gap_bed_file):
        a = line.split()
        if a[0] in valid_chrs:
            genome_bp -= int(a[2])-int(a[1])

    return genome_bp


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
