#!/usr/bin/env python
from optparse import OptionParser
import os, pdb, subprocess, sys, tempfile
import pysam
import bedtools, gff

################################################################################
# te_bam_instances.py
#
# Compute the number of aligned reads in a BAM file in individual instances of
# transposable elements.
#
# I'm thinking that provding a restricted GFF to the -t flag is the best
# way to focus on individual families or repeats.
#
# Bugs:
# -To intersect the BAM file with the repeats GFF, I need to use the -split
#   option to avoid intron intersections. However, then it splits each spliced
#   read so the same read could be counted twice if it intersects at both
#   junctions. I doubt this happens enough to worry about..
# -To filter the BAM file with a GFF, I require the read be fully within
#   the feature, but I need to relax this requirement for spliced reads
#   because otherwise their entire span is considered.
#
# WARNINGs:
# -Just don't use it for spliced transcripts as "-g", unless I have a control.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam_file>'
    parser = OptionParser(usage)    
    parser.add_option('-g', dest='filter_gff', help='Filter the TEs by overlap with genes in the given gff file [Default: %default]')
    parser.add_option('-s', dest='strand_split', default=False, action='store_true', help='Split statistics by strand [Default: %default]')
    parser.add_option('-t', dest='te_gff', default='%s/hg19.fa.out.tp.gff' % os.environ['MASK'])

    #parser.add_option('-f', dest='family', help='Limit to this family')
    #parser.add_option('-r', dest='repeat', help='Limit to this repeat')

    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide BAM file.')
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
        subprocess.call('intersectBed -a %s -b %s > %s' % (options.te_gff, filter_merged_bed_file, te_gff_file), shell=True)
        options.te_gff = te_gff_file

        # filter BAM
        bam_gff_fd, bam_gff_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
        bedtools.abam_f1(bam_files[i], filter_merged_bed_file, bam_gff_file)
        bam_file = bam_gff_file

        os.close(filter_merged_bed_fd)
        os.remove(filter_merged_bed_file)

    ############################################
    # count TE fragments
    ############################################
    fragments, te_fragments = count_te_fragments(bam_file, options.te_gff, options.strand_split)

    ############################################
    # print table
    ############################################
    for line in open(options.te_gff):
        a = line.split('\t')

        te_chrom = a[0]
        te_start = int(a[3])

        if options.strand_split:
            te_count = te_fragments.get((te_chrom, te_start, '+'), 0)
            te_pct = te_count / float(fragments)

            cols = (te_chrom, te_start, te_count, te_pct)
            print '%-5s  %9d  +  %6d  %9.2e' % cols

            te_count = te_fragments.get((te_chrom, te_start, '-'), 0)
            te_pct = te_count / float(fragments)

            cols = (te_chrom, te_start, te_count, te_pct)
            print '%-5s  %9d  -  %6d  %9.2e' % cols

        else:
            te_count = te_fragments.get((te_chrom, te_start), 0)
            te_pct = te_count / float(fragments)

            cols = (te_chrom, te_start, te_count, te_pct)
            print '%-5s  %9d  %6d  %9.2e' % cols

    ############################################
    # clean
    ############################################
    if options.filter_gff:
        os.close(te_gff_fd)
        os.remove(te_gff_file)

        os.close(bam_gff_fd)
        os.remove(bam_gff_file)


################################################################################
# count_te_fragments
#
# Count the number of fragments aligned to each individual TE.
################################################################################
def count_te_fragments(bam_file, te_gff, strand_split):
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

        if is_paired:
            read_inc = 0.5/multi_maps.get(a[3],1.0)
        else:
            read_inc = 1.0/multi_maps.get(a[3],1.0)

        te_chrom = a[12]
        te_start = int(a[15])
        te_kv = gff.gtf_kv(a[20])

        if strand_split:
            rstrand = a[5]
            tstrand = a[18]
            if rstrand == tstrand:
                te_key = (te_chrom, te_start, '+')
            else:
                te_key = (te_chrom, te_start, '-')
        else:
            te_key = (te_chrom, te_start)

        te_fragments[te_key] = te_fragments.get(te_key,0.0) + read_inc

    proc.communicate()

    return num_fragments, te_fragments


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
