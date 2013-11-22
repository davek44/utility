#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess, tempfile
import pysam

################################################################################
# count_reads.py
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
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide BAM file')

    print count(args[0])


################################################################################
# count
################################################################################
def count(bam_file, mapq=False):
    # initialize
    bam_count = 0.0

    # process
    for aligned_read in pysam.Samfile(bam_file, 'rb'):
        if mapq == False or aligned_read.mapq > 0:
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
# count_gtf
#
# WARNING:
#  I'm half assing by not considering strand information.
################################################################################
def count_gtf(bam_file, ref_gtf, mapq=False):
    ref_bam_fd, ref_bam_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])

    # intersect
    subprocess.call('intersectBed -abam %s -b %s > %s' % (bam_file, ref_gtf, ref_bam_file), shell=True)

    # count
    bam_count = count(ref_bam_file, mapq)

    # clean
    os.close(ref_bam_fd)
    os.remove(ref_bam_file)

    return bam_count


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
