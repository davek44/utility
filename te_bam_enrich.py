#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import binom
import gzip, os, random, subprocess, sys, tempfile
import pysam
import gff, stats

################################################################################
# te_bam_enrich.py
#
# Compute the enrichment of aligned reads in a BAM file in transposable
# element families.
#
# Slight bug:
# To intersect the BAM file with the repeats GFF, I need to use the -split
# option to avoid intron intersections. However, then it splits each spliced
# read so the same read could be counted twice if it intersects at both
# junctions. I doubt this happens much.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam file>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='control_bam_file', help='Control BAM file to paramterize null distribution [Default: %default]')
    parser.add_option('-g', dest='gff_file', help='Filter the TEs by overlap with genes in the given gff file [Default: %default]')
    parser.add_option('-m', dest='mapq', default=False, action='store_true', help='Consider only reads with mapq>0 [Default: %default]')
    parser.add_option('-r', dest='repeats_gff', default='%s/research/common/data/genomes/hg19/annotation/repeatmasker/hg19.fa.out.tp.gff' % os.environ['HOME'])
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide a gff file for the feature of interest.')
    else:
        bam_file = args[0]

    ############################################
    # GFF filter
    ############################################
    # filter TEs and read alignments by gff file
    if options.gff_file:
        # filter TE GFF
        te_gff_fd, te_gff_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
        subprocess.call('intersectBed -a %s -b %s > %s' % (options.repeats_gff, options.gff_file, te_gff_file), shell=True)
        options.repeats_gff = te_gff_file

        # filter BAM
        bam_gff_fd, bam_gff_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
        subprocess.call('intersectBed -abam %s -b %s > %s' % (bam_file, options.gff_file, bam_gff_file), shell=True)
        bam_file = bam_gff_file

        # filter control BAM
        if options.control_bam_file:
            cbam_gff_fd, cbam_gff_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
            subprocess.call('intersectBed -abam %s -b %s > %s' % (options.control_bam_file, options.gff_file, cbam_gff_file), shell=True)
            options.control_bam_file = cbam_gff_file

    ############################################
    # mapq filter
    ############################################
    if options.mapq:
        # filter BAM file for mapping quality
        bam_mapq_fd, bam_mapq_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
        bam_in = pysam.Samfile(bam_file, 'rb')
        bam_mapq_out = pysam.Samfile(bam_mapq_file, 'wb', template=bam_in)
        for aligned_read in bam_in:
            if aligned_read.mapq > 0:
                bam_mapq_out.write(aligned_read)
        bam_mapq_out.close()
        bam_file = bam_mapq_file

        # filter control BAM for mapping quality
        if options.control_bam_file:
            cbam_mapq_fd, cbam_mapq_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
            cbam_in = pysam.Samfile(options.control_bam_file, 'rb')
            cbam_mapq_out = pysam.Samfile(cbam_mapq_file, 'wb', template=cbam_in)
            for aligned_read in cbam_in:
                if aligned_read.mapq > 0:
                    cbam_mapq_out.write(aligned_read)
            cbam_mapq_out.close()
            options.control_bam_file = cbam_mapq_file

    ############################################
    # lengths
    ############################################
    # compute size of search space
    if options.gff_file:
        genome_length = count_gff(options.gff_file)
    else:
        genome_length = count_hg19()

    # estimate read length
    read_len = estimate_read_length(bam_file)

    # hash counted repeat genomic bp
    #te_lengths = measure_te(options.repeats_gff)
    te_lengths = te_target_size(options.repeats_gff, read_len)

    ############################################
    # count TE fragments
    ############################################
    fragments, te_fragments = count_te_fragments(bam_file, options.repeats_gff)
    if options.control_bam_file:
        control_fragments, control_te_fragments = count_te_fragments(options.control_bam_file, options.repeats_gff)

    ############################################
    # compute stats, print table
    ############################################
    for (rep,fam) in te_fragments:
        # parameterize null model
        if options.control_bam_file:
            te_p = control_te_fragments.get((rep,fam),1.0) / control_fragments
        else:
            te_p = float(te_lengths[(rep,fam)]) / genome_length

        # compute p-value of enrichment/depletion
        if te_fragments[(rep,fam)] > te_p*fragments:
            p_val = binom.sf(int(te_fragments[(rep,fam)])-1, int(fragments), te_p)
        else:
            p_val = binom.cdf(int(te_fragments[(rep,fam)]), int(fragments), te_p)

        # compute fold change
        if te_p*fragments > 0:
            fold = te_fragments[(rep,fam)]/(te_p*fragments)
        else:
            fold = 0

        cols = (rep, fam, te_lengths[(rep,fam)], te_fragments[(rep,fam)], te_p, fold, p_val)
        print '%-18s %-18s %10d %10d %10.2e %10.3f %10.2e' % cols

    ############################################
    # clean
    ############################################
    if options.gff_file:
        os.close(te_gff_fd)
        os.remove(te_gff_file)
        os.close(bam_gff_fd)
        os.remove(bam_gff_file)
        if options.control_bam_file:
            os.close(cbam_gff_fd)
            os.remove(cbam_gff_file)
    if options.mapq:
        os.close(bam_mapq_fd)
        os.remove(bam_mapq_file)
        if options.control_bam_file:
            os.close(cbam_mapq_fd)
            os.remove(cbam_mapq_file)


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
# count_te_fragments
#
# Count the number of fragments aligned to each TE family.
################################################################################
def count_te_fragments(bam_file, te_gff):
    # count fragments and hash multi-mappers
    num_fragments = 0
    multi_maps = {}
    paired_poll = {False:0, True:0}
    for aligned_read in pysam.Samfile(bam_file, 'rb'):
        if aligned_read.is_paired:
            num_fragments += 0.5/aligned_read.opt('NH')
        else:
            num_fragments += 1.0/aligned_read.opt('NH')

        if aligned_read.opt('NH') > 1:
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
    te_fragments = {}
    proc = subprocess.Popen('intersectBed -split -wo -bed -abam %s -b %s' % (bam_file, te_gff), shell=True, stdout=subprocess.PIPE)
    for line in proc.stdout:
        a = line.split('\t')
        te_kv = gff.gtf_kv(a[14])

        if is_paired:
            read_inc = 0.5/multi_maps.get(a[3],1.0)
        else:
            read_inc = 1.0/multi_maps.get(a[3],1.0)

        te_fragments[(te_kv['repeat'],te_kv['family'])] = te_fragments.get((te_kv['repeat'],te_kv['family']),0.0) + read_inc
        te_fragments[('*',te_kv['family'])] = te_fragments.get(('*',te_kv['family']),0.0) + read_inc
        te_fragments[('*','*')] = te_fragments.get(('*','*'),0.0) + read_inc

    proc.communicate()

    return num_fragments, te_fragments


################################################################################
# estimate_read_length
#
# Compute mean read length by sampling the first N reads.
################################################################################
def estimate_read_length(bam_file):
    samples = 10000
    s = 0
    read_lengths = []
    for aligned_read in pysam.Samfile(bam_file, 'rb'):
        read_lengths.append(aligned_read.rlen)
        s += 1
        if s >= samples:
            break
    return int(0.5+stats.mean(read_lengths))
        

################################################################################
# measure_te
#
# Hash the number of bp covered by various repeats in the RepeatMasker gff file
# and the lincRNA gtf file.
#
# Note:
#  -This version is obsolete because it doesn't consider that overlapping
#   reads need not be inside the TE and the true range depends on the read
#   length. It's been replaced by te_target_size.
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
# te_target_size
#
# Measure the overlap target area for each TE for the given read length.
################################################################################
def te_target_size(te_gff, read_len):
    te_bp = {}
    active_te_intervals = {}

    for line in open(te_gff):
        a = line.split('\t')

        kv = gff.gtf_kv(a[8])
        rep = kv['repeat']
        fam = kv['family']

        chrom = a[0]
        start = int(a[3])
        end = int(a[4])

        # process closed intervals
        for arep, afam in active_te_intervals.keys():
            achrom, astart, aend = active_te_intervals[(arep,afam)]

            if achrom != chrom or aend + read_len < start:
                # add
                te_bp[(arep,afam)] = te_bp.get((arep,afam),0) + aend - astart + 1 + read_len
                # close
                del active_te_intervals[(arep,afam)]

        # update/add te
        if (rep,fam) in active_te_intervals:
            achrom, astart, aend = active_te_intervals[(rep,fam)]
            active_te_intervals[(rep,fam)] = (chrom, min(astart,start), max(aend, end))
        else:
            active_te_intervals[(rep,fam)] = (chrom, start, end)

        if ('*',fam) in active_te_intervals:
            achrom, astart, aend = active_te_intervals[('*',fam)]
            active_te_intervals[('*',fam)] = (chrom, min(astart,start), max(aend, end))
        else:
            active_te_intervals[('*',fam)] = (chrom, start, end)

        if ('*','*') in active_te_intervals:
            achrom, astart, aend = active_te_intervals[('*','*')]
            active_te_intervals[('*','*')] = (chrom, min(astart,start), max(aend, end))
        else:
            active_te_intervals[('*','*')] = (chrom, start, end)

    # close remaining
    for arep, afam in active_te_intervals.keys():
        achrom, astart, aend = active_te_intervals[(arep,afam)]

        # add
        te_bp[(arep,afam)] = te_bp.get((arep,afam),0) + aend - astart + 1 + read_len

    return te_bp


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
