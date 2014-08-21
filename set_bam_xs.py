#!/usr/bin/env python
from optparse import OptionParser
import os
import pysam

################################################################################
# set_bam_xs.py
#
# Set the XS tag properly in a BAM file. Currently assumes first-strand because
# that's all that I've seen so far.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam>'
    parser = OptionParser(usage)
    parser.add_option('-a', dest='alter_bam', help='Print alterations to this file [Default: %default]')
    parser.add_option('-o', dest='bam_out_file', help='Output BAM file')
    parser.add_option('-r', dest='reverse_protocol', default=False, action='store_true', help='For fr-secondstrand [Default: %default]')
    parser.add_option('-s', dest='strip', action='store_true', default=False, help='Strip XS tags for unspliced reads [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide BAM file')
    else:
        bam_in_file = args[0]

    if options.bam_out_file == None:
        b, e = os.path.splitext(bam_in_file)
        options.bam_out_file = b+'_xs'+e

    bam_in = pysam.Samfile(bam_in_file, 'rb')
    bam_out = pysam.Samfile(options.bam_out_file, 'wb', template=bam_in)

    if options.alter_bam:
        alter_bam_out = pysam.Samfile(options.alter_bam, 'wb', template=bam_in)

    for aligned_read in bam_in:
        # just strip the tag
        if options.strip:
            # ... of unspliced reads
            if not spliced(aligned_read):
                # ... with a tag
                try:
                    xs_tag = aligned_read.opt('XS')
                except:
                    xs_tag = None

                if xs_tag:
                    # remove tag
                    rm_xs(aligned_read)

                    # fix CP tag
                    fix_cp(aligned_read)

            # output
            bam_out.write(aligned_read)

        # set the tag properly
        else:
            # determine XS
            if aligned_read.is_paired:
                if aligned_read.is_read1:
                    if aligned_read.is_reverse:
                        new_xs = '+'
                    else:
                        new_xs = '-'
                else:
                    if aligned_read.is_reverse:
                        new_xs = '-'
                    else:
                        new_xs = '+'
            else:
                if aligned_read.is_reverse:
                    new_xs = '-'
                else:
                    new_xs = '+'

            if options.reverse_protocol:
                if new_xs == '+':
                    new_xs = '-'
                else:
                    new_xs = '+'

            # toss it if the splicing strand differs
            if not splice_disagree(aligned_read, new_xs):
                # get old
                try:
                    old_xs = aligned_read.opt('XS')
                except:
                    old_xs = None

                if old_xs == new_xs:
                    # write as is
                    bam_out.write(aligned_read)

                else:
                    # remove existing XS tag
                    rm_xs(aligned_read)

                    # set XS tag
                    aligned_read.tags = aligned_read.tags + [('XS',new_xs)]

                    # fix CP tag
                    fix_cp(aligned_read)

                    # output
                    bam_out.write(aligned_read)

                    if options.alter_bam:
                        alter_bam_out.write(aligned_read)

    bam_in.close()
    bam_out.close()
    if options.alter_bam:
        alter_bam_out.close()
    

################################################################################
# fix_cp
#
# The CP tag is read as a float rather than int, and then improperly handled
# by the samtools library.
################################################################################
def fix_cp(aligned_read):
    cp_i = 0
    while cp_i < len(aligned_read.tags) and aligned_read.tags[cp_i][0] != 'CP':
        cp_i += 1

    if cp_i < len(aligned_read.tags):
        cp_int = int(aligned_read.tags[cp_i][1])
        aligned_read.tags = aligned_read.tags[:cp_i] + aligned_read.tags[cp_i+1:] + [('CP',cp_int)]


################################################################################
# rm_xs
#
# Remove the XS tag from the AlignedRead object
################################################################################
def rm_xs(aligned_read):
    xs_i = 0
    while xs_i < len(aligned_read.tags) and aligned_read.tags[xs_i][0] != 'XS':
        xs_i += 1

    if xs_i < len(aligned_read.tags):
        aligned_read.tags = aligned_read.tags[:xs_i] + aligned_read.tags[xs_i+1:]


################################################################################
# splice_disagree
#
# Return true if the read is spliced and the current XS tag disagrees with the
# new one.
################################################################################
def splice_disagree(aligned_read, new_xs):
    if not spliced(aligned_read):
        return False
    else:
        return new_xs != aligned_read.opt('XS')


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
