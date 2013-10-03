#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import binom
import gzip, os, pdb, random, subprocess, sys, tempfile
import pysam
import bedtools, gff, stats

################################################################################
# te_bam_enrich.py
#
# Compute the enrichment of aligned reads in a BAM file in transposable
# element families.
#
# Bugs:
# -To intersect the BAM file with the repeats GFF, I need to use the -split
#   option to avoid intron intersections. However, then it splits each spliced
#   read so the same read could be counted twice if it intersects at both
#   junctions. I doubt this happens much.
# -To filter the BAM file with a GFF, I require the read be fully within
#   the feature, but I need to relax this requirement for spliced reads
#   because otherwise their entire span is considered.
# -Spliced transcripts as "-g" will be undercounted by "count_bed" because
#   I assume each region is a unit which the read must be contained within.
#
# WARNINGs:
# -Just don't use it for spliced transcripts as "-g", unless I have a control.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam file>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='control_bam_file', help='Control BAM file to paramterize null distribution [Default: %default]')
    parser.add_option('-g', dest='filter_gff', help='Filter the TEs by overlap with genes in the given gff file [Default: %default]')
    parser.add_option('-m', dest='mapq', default=False, action='store_true', help='Consider only reads with mapq>0 [Default: %default]')
    parser.add_option('-r', dest='repeats_gff', default='%s/hg19.fa.out.tp.gff' % os.environ['MASK'])
    parser.add_option('-s', dest='strand_split', default=False, action='store_true', help='Split statistics by strand [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide a BAM file.')
    else:
        bam_file = args[0]

    ############################################
    # GFF filter
    ############################################
    # filter TEs and read alignments by gff file
    if options.filter_gff:
        filter_merged_bed_fd, filter_merged_bed_file = tempfile.mkstemp()
        subprocess.call('sortBed -i %s | mergeBed -i - > %s' % (options.filter_gff, filter_merged_bed_file), shell=True)

        # filter TE GFF
        te_gff_fd, te_gff_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
        subprocess.call('intersectBed -a %s -b %s > %s' % (options.repeats_gff, filter_merged_bed_file, te_gff_file), shell=True)
        options.repeats_gff = te_gff_file

        # filter BAM
        bam_gff_fd, bam_gff_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
        bedtools.abam_f1(bam_file, filter_merged_bed_file, bam_gff_file)
        bam_file = bam_gff_file

        # filter control BAM
        if options.control_bam_file:
            cbam_gff_fd, cbam_gff_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
            bedtools.abam_f1(options.control_bam_file, filter_merged_bed_file, cbam_gff_file)
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
    # estimate read length
    read_len = estimate_read_length(bam_file)

    # compute size of search space
    if options.filter_gff:
        genome_length = count_bed(filter_merged_bed_file, read_len)
    else:
        genome_length = count_hg19()

    # hash counted repeat genomic bp
    #te_lengths = measure_te(options.repeats_gff)
    if options.filter_gff:
        te_lengths = te_target_size_bed(options.repeats_gff, filter_merged_bed_file, read_len)
    else:
        te_lengths = te_target_size(options.repeats_gff, read_len)

    ############################################
    # count TE fragments
    ############################################
    fragments, te_fragments = count_te_fragments(bam_file, options.repeats_gff, options.strand_split)
    all_tes = set(te_fragments.keys())

    if options.control_bam_file:
        control_fragments, control_te_fragments = count_te_fragments(options.control_bam_file, options.repeats_gff, options.strand_split)

        # add control pseudocounts
        all_tes |= control_te_fragments.keys())
        for (rep,fam) in all_tes:
            if rep != '*':
                control_te_fragments[(rep,fam)] = control_te_fragments.get((rep,fam),0) + 1
                control_te_fragments[('*',fam)] = control_te_fragments.get(('*',fam),0) + 1
                control_te_fragments[('*','*')] = control_te_fragments.get(('*','*'),0) + 1
                control_fragments += 1
            if not (rep,fam) in te_fragments:
                te_fragments[(rep,fam)] = 0

    ############################################
    # compute stats, print table
    ############################################
    for (rep,fam) in all_tes:
        # parameterize null model
        if options.control_bam_file:
            te_p = control_te_fragments[(rep,fam)] / control_fragments
        else:
            if options.strand_split:
                te_p = float(te_lengths[(rep[:-1],fam)]) / (2*genome_length)
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

        if options.strand_split:
            te_len = te_lengths[(rep[:-1],fam)]
        else:
            te_len = te_lengths[(rep,fam)]
            
        cols = (rep, fam, te_len, te_fragments[(rep,fam)], te_p*fragments, fold, p_val)
        print '%-18s %-18s %10d %10.1f %10.1f %10.3f %10.2e' % cols

    ############################################
    # clean
    ############################################
    if options.filter_gff:
        os.close(filter_merged_bed_fd)
        os.remove(filter_merged_bed_file)
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
# count_bed
#
# Count the number of bp in the filter merged BED file.
################################################################################
def count_bed(bed_file, read_len):
    bp = 0    
    for line in open(bed_file):
        a = line.split()
        bp += int(a[2]) - int(a[1]) - read_len + 1
    return bp


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
def count_te_fragments(bam_file, te_gff, strand_split=False):
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
        te_kv = gff.gtf_kv(a[20])

        rep = te_kv['repeat']
        fam = te_kv['family']

        if is_paired:
            read_inc = 0.5/multi_maps.get(a[3],1.0)
        else:
            read_inc = 1.0/multi_maps.get(a[3],1.0)

        rep_star = '*'
        if strand_split:
            rstrand = a[5]
            tstrand = a[18]
            if rstrand == tstrand:
                rep += '+'
                rep_star += '+'
            else:
                rep += '-'
                rep_star += '-'

        te_fragments[(rep,fam)] = te_fragments.get((rep,fam),0.0) + read_inc
        te_fragments[(rep_star,fam)] = te_fragments.get((rep_star,fam),0.0) + read_inc
        te_fragments[(rep_star,'*')] = te_fragments.get((rep_star,'*'),0.0) + read_inc

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
#  -This version is no longer used to parameterize the null model  because it
#   doesn't consider that overlapping reads need not be inside the TE and the
#   true range depends on the read length. It's been replaced by te_target_size.
#   However, it's still used by a few different places.
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
# te_target_size_bed
#
# Measure the overlap target area for each TE within each BED interval for the
# given read length.
################################################################################
def te_target_size_bed(te_gff, ref_bed, read_len):
    # hash TE intervals by BED region
    bed_te_intervals = {}
    p = subprocess.Popen('intersectBed -wo -a %s -b %s' % (ref_bed, te_gff), shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        a = line.split('\t')

        bchrom = a[0]
        bstart = int(a[1])
        bend = int(a[2])
        bid = (bchrom,bstart)

        rep_kv = gff.gtf_kv(a[11])
        rep = rep_kv['repeat']
        fam = rep_kv['family']

        tstart = int(a[6])
        tend = int(a[7])

        ostart = max(bstart, tstart)
        oend = min(bend, tend)

        if not bid in bed_te_intervals:
            bed_te_intervals[bid] = {}
        bed_te_intervals[bid].setdefault((rep,fam),[]).append((ostart,oend))
        bed_te_intervals[bid].setdefault(('*',fam),[]).append((ostart,oend))
        bed_te_intervals[bid].setdefault(('*','*'),[]).append((ostart,oend))        

    p.communicate()

    target_size = {}
    for bid in bed_te_intervals:
        bchrom, bstart = bid        

        for te in bed_te_intervals[bid]:
            bt_intervals = bed_te_intervals[bid][te]
            bt_intervals.sort()

            # merge intervals, limited at the start by the BED region's start
            merged_intervals = [(max(bstart, bt_intervals[0][0]-read_len+1), bt_intervals[0][1])]
            for i in range(1,len(bt_intervals)):
                start1, end1 = merged_intervals[-1]
                start2, end2 = bt_intervals[i]

                if end1+1 < start2-read_len+1:
                    merged_intervals.append((start2-read_len+1,end2))
                else:
                    merged_intervals[-1] = (start1, end2)

            # sum
            target_size[te] = target_size.get(te,0) + sum([e-s+1 for (s,e) in merged_intervals])

    return target_size


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
