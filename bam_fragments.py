#!/usr/bin/env python
from __future__ import print_function
from optparse import OptionParser
import os, subprocess, tempfile
import pysam

################################################################################
# bam_fragments.py
#
# Count the reads in a BAM file, by relying on the multimap count tag.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam>'
    parser = OptionParser(usage)
    #parser.add_option('-d', dest='debug', default=False, action='store_true', help='Print debug output to read_counts.txt')
    parser.add_option('-g', dest='gff_file', help='Count fragments overlapping the GFF file')
    parser.add_option('-m', dest='filter_mapq', default=False, action='store_true', help='Filter out mapq=0 alignments [Default: %default]')
    parser.add_option('-t', dest='gtf_file', help='Count fragments overlapping the GTF file')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide BAM file')
    else:
        bam_file = args[0]

    if options.gff_file:
        print(count_gff(bam_file, options.gff_file, options.filter_mapq))
    elif options.gtf_file:
        print(count_gtf(bam_file, options.gtf_file, options.filter_mapq))
    else:
        print(count(bam_file))


################################################################################
# count
#
# Note: segemehl writes NH as the total possible alignments.
################################################################################
def count(bam_file, filter_mapq=False):
    # initialize
    bam_count = 0.0

    # process
    bam_in = pysam.Samfile(bam_file, 'rb')
    for aligned_read in bam_in:
        if filter_mapq == False or aligned_read.mapq > 0:
            nh_tag = 1
            if aligned_read.has_tag('NH'):
                nh_tag = aligned_read.get_tag('NH')

            if aligned_read.is_paired:
                bam_count += 0.5/nh_tag
            else:
                bam_count += 1.0/nh_tag

    return bam_count


################################################################################
# count_gff
#
# Count reads in the bam_file that intersect the GFF file.
#
# Note: I need to use -split -bed because otherwise spliced reads overlap
#       gff entries in between their aligned segments.
################################################################################
def count_gff(bam_file, gff_file, filter_mapq=False):
    bam_in = pysam.Samfile(bam_file, 'rb')

    # segemehl
    aligner = bam_in.header['PG'][0]['ID']
    if aligner[-8:] == 'segemehl':
        multimap_max = segemehl_multimap_max(bam_in)

    # hash multi-mappers
    multi_maps = {}
    paired_poll = {False:0, True:0}
    for aligned_read in bam_in:
        if aligned_read.get_tag('NH') > 1:
            nh_tag = aligned_read.get_tag('NH')
            if aligner == 'segemehl':
                multi_maps[aligned_read.qname] = min(nh_tag, multimap_max)
            else:
                multi_maps[aligned_read.qname] = nh_tag

        paired_poll[aligned_read.is_paired] += 1

    # guess paired-ness
    if paired_poll[True] > 0 and paired_poll[False] > 0:
        print('Paired-ness of the reads is ambiguous', file=sys.stderr)
    if paired_poll[True] > paired_poll[False]:
        is_paired = True
    else:
        is_paired = False

    # intersect and count
    bam_count = 0.0
    p = subprocess.Popen('bedtools intersect -split -u -bed -abam %s -b %s' % (bam_file,gff_file), shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        a = line.split('\t')
        if is_paired:
            bam_count += 0.5/multi_maps.get(a[3],1.0)
        else:
            bam_count += 1.0/multi_maps.get(a[3],1.0)
    p.communicate()

    return bam_count


################################################################################
# count_gtf
#
# WARNING:
#  -I'm half assing by not considering strand information.
#  -Also, I'm skipping the -split -bed stuff because spliced reads overlapping
#    genes are extemely likely to count anyway.
################################################################################
def count_gtf(bam_file, ref_gtf, filter_mapq=False):
    ref_bam_fd, ref_bam_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])

    # intersect
    subprocess.call('bedtools intersect -abam %s -b %s > %s' % (bam_file, ref_gtf, ref_bam_file), shell=True)

    # count
    bam_count = count(ref_bam_file, filter_mapq)

    # clean
    os.close(ref_bam_fd)
    os.remove(ref_bam_file)

    return bam_count


################################################################################
# segemehl_multimap_max
#
# Grab the -r multimap max parameter from a segemehl alignment command.
################################################################################
def segemehl_multimap_max(bam_in):
    align_cmd = bam_in.header['PG'][0]['CL']
    optr_i = align_cmd.find('-r')

    if optr_i == -1:
        # what's the default?
        multimap_max = 20
    else:
        optr_start = optr_i+2

        # get past possible initial spaces
        optr_end = optr_start
        while align_cmd[optr_end] == ' ':
            optr_end += 1

        # find the end of the number
        while align_cmd[optr_end].isdigit():
            optr_end += 1

        # grab the number
        multimap_max = int(align_cmd[optr_start:optr_end])

    return multimap_max


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
