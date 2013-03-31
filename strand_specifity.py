#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess
import pysam

################################################################################
# strand_specificity.py
#
# Print information relevant to determining the strand specificity of the
# sequencing in a BAM file using a TopHat junctions.bed file.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam> <junctions.bed>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide BAM file and junctions.bed file.')
    else:
        bam_file = args[0]
        juncs_bed_file = args[1]

    # filter junctions for forward only
    subprocess.call('awk \'$6 == "+"\' %s > junctions_fwd.bed' % juncs_bed_file, shell=True)

    # intersect BAM with forward junctions
    subprocess.call('intersectBed -s -abam %s -b junctions_fwd.bed > fwd.bam' % bam_file, shell=True)

    # count first/second reads
    first = 0
    second = 0
    for aligned_read in pysam.Samfile('fwd.bam'):
        if aligned_read.is_proper_pair:
            if aligned_read.is_read1:
                first += 1
            else:
                second += 1

    print 'Read1\'s aligning + and intersecting + junctions: %9d' % first
    print 'Read2\'s aligning + and intersecting + junctions: %9d' % second

    os.remove('junctions_fwd.bed')
    os.remove('fwd.bam')


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
