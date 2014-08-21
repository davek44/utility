#!/usr/bin/env python
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
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide BAM file')
    else:
        bam_file = args[0]

    if options.gff_file:
        print count_gff(bam_file, options.gff_file)
    else:
        print count(bam_file)


################################################################################
# count
################################################################################
def count(bam_file, filter_mapq=False):
    # initialize
    bam_count = 0.0

    # process
    for aligned_read in pysam.Samfile(bam_file, 'rb'):
        if filter_mapq == False or aligned_read.mapq > 0:
            try:
                nh_tag = aligned_read.opt('NH')
            except:
                nh_tag = 1

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
    # hash multi-mappers
    multi_maps = {}
    paired_poll = {False:0, True:0}
    for aligned_read in pysam.Samfile(bam_file, 'rb'):
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

    # intersect and count
    bam_count = 0.0
    p = subprocess.Popen('intersectBed -split -wo -bed -abam %s -b %s' % (bam_file,gff_file), shell=True, stdout=subprocess.PIPE)
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
    subprocess.call('intersectBed -abam %s -b %s > %s' % (bam_file, ref_gtf, ref_bam_file), shell=True)

    # count
    bam_count = count(ref_bam_file, filter_mapq)

    # clean
    os.close(ref_bam_fd)
    os.remove(ref_bam_file)

    return bam_count


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
