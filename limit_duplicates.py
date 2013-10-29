#!/usr/bin/env python
from optparse import OptionParser
import pdb, os, random, sys, subprocess, tempfile
import pysam

################################################################################
# limit_duplicates.py
#
# Accept a sorted BAM file of unpaired reads as input and max out the number of
# reads that can occur at one single position.
#
# For now, I'm just randomly choosing a read rather than doing any work to
# use the "best" read by some definition or make the choice consistent across
# multiple mapping reads. 
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam file>'
    parser = OptionParser(usage)
    parser.add_option('-d', dest='dup_t', type='int', default=10, help='Number of duplicates at which removal begins [Default: %default]')
    parser.add_option('-o', dest='output_bam', default='out.bam', help='Ouput BAM file')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide bam file')
    else:
        bam_file = args[0]

    bam_in = pysam.Samfile(bam_file, 'rb')

    bam_out_fd, bam_out_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
    bam_out = pysam.Samfile(bam_out_file, 'wb', header=bam_in.header)

    forward_position_reads = []
    forward_position = None
    reverse_position_read_hash = {}

    for aligned_read in bam_in:

        if aligned_read.is_reverse:
            ar_pos = (aligned_read.tid, aligned_read.aend)

            # add this read
            reverse_position_read_hash.setdefault(ar_pos,[]).append(aligned_read)

            # consider closing previous positions
            for reverse_position in reverse_position_read_hash.keys():
                if reverse_position[0] != aligned_read.tid or reverse_position[1] < aligned_read.pos:
                    # output
                    sample_write(bam_out, reverse_position_read_hash[reverse_position], options.dup_t)

                    # close
                    del reverse_position_read_hash[reverse_position]                    
            
        else:
            ar_pos = (aligned_read.tid, aligned_read.pos)

            if forward_position == None or forward_position == ar_pos:
                # same position, keep adding
                forward_position_reads.append(aligned_read)
            else:
                # new position, output
                sample_write(bam_out, forward_position_reads, options.dup_t)

                # close and update
                forward_position_reads = [aligned_read]

            forward_position = ar_pos

    # print remaining
    sample_write(bam_out, forward_position_reads, options.dup_t)
    for reverse_position in reverse_position_read_hash:
        sample_write(bam_out, reverse_position_read_hash[reverse_position], options.dup_t)

    bam_in.close()
    bam_out.close()

    # re-sort
    output_bam_prefix = os.path.splitext(options.output_bam)[0]
    subprocess.call('samtools sort %s %s' % (bam_out_file, output_bam_prefix), shell=True)

    os.close(bam_out_fd)
    os.remove(bam_out_file)


################################################################################
# sample_write
#
# Sample dup_t reads from the list and output.
################################################################################
def sample_write(bam_out, position_reads, dup_t):
    if len(position_reads) > dup_t:
        position_reads = random.sample(position_reads, dup_t)
    for sampled_aligned_read in position_reads:
        bam_out.write(sampled_aligned_read)

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
