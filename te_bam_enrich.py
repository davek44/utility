#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import binom
import glob, os, pdb, random, subprocess, sys, tempfile
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
    usage = 'usage: %prog [options] <bam_file,bam_file2,...>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='control_bam_files', help='Control BAM file to paramterize null distribution [Default: %default]')
    parser.add_option('-f', dest='filter_gff', help='Filter the TEs by overlap with genes in the given gff file [Default: %default]')
    parser.add_option('-g', dest='genome', default='HG19', help='Genome directory to obtain lengths from [Default: %default]')
    parser.add_option('-m', dest='mapq', default=False, action='store_true', help='Consider only reads with mapq>0 [Default: %default]')
    parser.add_option('-r', dest='repeats_gff', default='%s/hg19.fa.out.tp.gff' % os.environ['MASK'])
    parser.add_option('-s', dest='strand_split', default=False, action='store_true', help='Split statistics by strand [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide BAM files.')
    else:
        bam_files = args[0].split(',')

    control_bam_files = []
    if options.control_bam_files:
        control_bam_files = options.control_bam_files.split(',')

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
        bam_gff_fds = [None]*len(bam_files)
        bam_gff_files = [None]*len(bam_files)
        for i in range(len(bam_files)):
            bam_gff_fds[i], bam_gff_files[i] = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
            bedtools.abam_f1(bam_files[i], filter_merged_bed_file, bam_gff_files[i])
            bam_files[i] = bam_gff_files[i]

        # filter control BAM
        if control_bam_files:
            cbam_gff_fds = [None]*len(control_bam_files)
            cbam_gff_files = [None]*len(control_bam_files)
            for i in range(len(control_bam_files)):
                cbam_gff_fds[i], cbam_gff_files[i] = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
                bedtools.abam_f1(control_bam_files[i], filter_merged_bed_file, cbam_gff_files[i])
                control_bam_files[i] = cbam_gff_files[i]

    ############################################
    # lengths
    ############################################
    # estimate read length (just averaging across replicates for now)
    read_lens = []
    for bam_file in bam_files:
        read_lens.append(estimate_read_length(bam_file))
    read_len = stats.mean(read_lens)

    # compute size of search space
    if options.filter_gff:
        genome_length = count_bed(filter_merged_bed_file, read_len)
    else:
        genome_length = count_genome(options.genome)

    # hash counted repeat genomic bp
    if options.filter_gff:
        te_lengths = te_target_size_bed(options.repeats_gff, filter_merged_bed_file, read_len)
    else:
        te_lengths = te_target_size(options.repeats_gff, read_len)

    ############################################
    # count TE fragments
    ############################################
    fragments = []
    te_fragments = []
    for bam_file in bam_files:
        rep_fragments, rep_te_fragments = count_te_fragments(bam_file, options.repeats_gff, options.strand_split)
        fragments.append(rep_fragments)
        te_fragments.append(rep_te_fragments)

    if control_bam_files:
        control_fragments = []
        control_te_fragments = []
        for control_bam_file in control_bam_files:
            rep_fragments, rep_te_fragments = count_te_fragments(control_bam_file, options.repeats_gff, options.strand_split)
            control_fragments.append(rep_fragments)
            control_te_fragments.append(rep_te_fragments)

    ############################################
    # combine replicates into fragment rates
    ############################################
    te_fragment_rates = {}
    for (rep,fam) in te_lengths:
        if options.strand_split:
            # positive
            rate_list = [te_fragments[i].get((rep+'+',fam),1)/float(fragments[i]) for i in range(len(bam_files))]
            te_fragment_rates[(rep+'+',fam)] = stats.geo_mean(rate_list)
            # negative
            rate_list = [te_fragments[i].get((rep+'-',fam),1)/float(fragments[i]) for i in range(len(bam_files))]
            te_fragment_rates[(rep+'-',fam)] = stats.geo_mean(rate_list)
        else:
            rate_list = [te_fragments[i].get((rep,fam),1)/float(fragments[i]) for i in range(len(bam_files))]
            te_fragment_rates[(rep,fam)] = stats.geo_mean(rate_list)

    if control_bam_files:
        control_te_fragment_rates = {}
        for te in te_fragment_rates:
            rate_list = [control_te_fragments[i].get(te,1)/float(control_fragments[i]) for i in range(len(control_bam_files))]
            control_te_fragment_rates[te] = stats.geo_mean(rate_list)

    ############################################
    # compute stats, print table
    ############################################
    for (rep,fam) in te_fragment_rates:
        # compute TE length
        if options.strand_split:
            te_len = te_lengths[(rep[:-1],fam)]
        else:
            te_len = te_lengths[(rep,fam)]

        # parameterize null model
        if options.control_bam_files:
            null_rate = control_te_fragment_rates[(rep,fam)]
        else:
            if options.strand_split:
                null_rate = float(te_lengths[(rep[:-1],fam)]) / (2*genome_length)
            else:
                null_rate = float(te_lengths[(rep,fam)]) / genome_length

        # compute fragment counts
        count = te_fragment_rates[(rep,fam)]*sum(fragments)
        null_count = null_rate*sum(fragments)

        # compute fold change
        if null_rate > 0:
            fold = te_fragment_rates[(rep,fam)]/null_rate
        else:
            fold = 0

        # compute p-value of enrichment/depletion
        p_val = 1.0
        for i in range(len(bam_files)):
            if te_fragment_rates[(rep,fam)] > null_rate:
                p_val *= binom.sf(int(te_fragments[i].get((rep,fam),1))-1, int(fragments[i]), null_rate)
            else:
                p_val *= binom.cdf(int(te_fragments[i].get((rep,fam),1)), int(fragments[i]), null_rate)

        cols = (rep, fam, te_len, count, null_count, fold, p_val)
        print '%-18s %-18s %10d %10.1f %10.1f %10.3f %10.2e' % cols

    ############################################
    # clean
    ############################################
    if options.filter_gff:
        os.close(filter_merged_bed_fd)
        os.remove(filter_merged_bed_file)
        os.close(te_gff_fd)
        os.remove(te_gff_file)

        for i in range(len(bam_files)):
            os.close(bam_gff_fds[i])
            os.remove(bam_gff_files[i])

        if options.control_bam_files:
            for i in range(len(control_bam_files)):
                os.close(cbam_gff_fds[i])
                os.remove(cbam_gff_files[i])


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
# count_genome
#
# Count the number of bp in the genome where TEs could be.
################################################################################
def count_genome(genome):
    chrom_sizes_file = glob.glob('%s/assembly/*.genome' % os.environ[genome])[0]

    valid_chrs = ['chr%d' % c for c in range(1,23)] + ['chrX','chrY']

    genome_bp = 0
    for line in open(chrom_sizes_file):
        a = line.split()
        if len(a) > 0:
            if a[0].find('random') != -1:
                continue
            if a[0].find('chrUn') != -1:
                continue
            if a[0].find('hap') != -1:
                continue

            genome_bp += int(a[1])

    gap_bed_files = glob.glob('%s/assembly/*gaps.bed' % os.environ[genome])
    if len(gap_bed_files):
        gap_bed_file = gap_bed_files[0]
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
        try:
            nh = aligned_read.get_tag('NH')
        except KeyError:
            nh = 1

        if aligned_read.is_paired:
            num_fragments += 0.5/nh
        else:
            num_fragments += 1.0/nh

        if nh > 1:
            multi_maps[aligned_read.qname] = nh

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
