#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess, tempfile
import pysam

################################################################################
# bedtools.py
#
#
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()
    

################################################################################
# abam_f1
#
# Intersect the BAM file with the BED file using the "-f 1" option, but correct
# for the loss of spliced reads.
################################################################################
def abam_f1(bam_file, bed_file, out_file):
    ############################################
    # divide BAM by splicing
    ############################################
    spliced_bam_fd, spliced_bam_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
    unspliced_bam_fd, unspliced_bam_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])

    # open BAMs
    bam_in = pysam.Samfile(bam_file, 'rb')
    spliced_bam_out = pysam.Samfile(spliced_bam_file, 'wb', template=bam_in)
    unspliced_bam_out = pysam.Samfile(unspliced_bam_file, 'wb', template=bam_in)

    # divide
    for aligned_read in bam_in:
        if spliced(aligned_read):
            spliced_bam_out.write(aligned_read)
        else:
            unspliced_bam_out.write(aligned_read)

    # close BAMs
    bam_in.close()
    spliced_bam_out.close()
    unspliced_bam_out.close()

    ############################################
    # intersect and merge
    ############################################
    spliced_is_bam_fd, spliced_is_bam_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
    unspliced_is_bam_fd, unspliced_is_bam_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])

    subprocess.call('intersectBed -f 1 -abam %s -b %s > %s' % (unspliced_bam_file, bed_file, unspliced_is_bam_file), shell=True)
    subprocess.call('intersectBed -abam %s -b %s > %s' % (spliced_bam_file, bed_file, spliced_is_bam_file), shell=True)

    subprocess.call('samtools merge -f %s %s %s' % (out_file, unspliced_is_bam_file, spliced_is_bam_file), shell=True)

    ############################################
    # clean
    ############################################
    os.close(spliced_bam_fd)
    os.remove(spliced_bam_file)
    os.close(unspliced_bam_fd)
    os.remove(unspliced_bam_file)
    os.close(spliced_is_bam_fd)
    os.remove(spliced_is_bam_file)
    os.close(unspliced_is_bam_fd)
    os.remove(unspliced_is_bam_file)
    

################################################################################
# spliced
#
# Return true if the read is spliced.
################################################################################
def spliced(aligned_read):
    spliced = False
    for code,size in aligned_read.cigar:
        if code == 3:
            spliced = True
    return spliced

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
