#!/usr/bin/env python
from optparse import OptionParser
import copy, gzip
import gff

################################################################################
# rm2gff_remove_overlaps.py
#
# Convert RepeatMasker .out format to gff to and trim overlapping annotations.
#
# I've found that the script doesn't totally work, especially for simple
# repeats. A better method would be to take each connected component of
# overlapping repeats and find the maximum weight non-overlapping set.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <rm out>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide RepeatMasker .out file')
    else:
        if args[0][-2:] == 'gz':
            rm_in = gzip.open(args[0])
        else:
            rm_in = open(args[0])

    for i in range(4):
        line = rm_in.readline()

    last_chrom = ''
    last_a = []
    while line:
        if last_a:
            last_score = int(last_a[0])
            last_chrom = last_a[4]
            last_start = int(last_a[5])
            last_end = int(last_a[6])

        this_a = line.split()
        this_score = int(this_a[0])
        this_chrom = this_a[4]
        this_start = int(this_a[5])
        this_end = int(this_a[6])
        
        if last_chrom != this_chrom or this_start > last_end:
            # print the prior line
            if last_a:
                print gff_line(last_a)

            # move this to last
            last_a = this_a

        else:
            if this_start < last_start:
                # nasty 3 way overlap, just move this
                this_a[5] = str(last_start)
                this_start = last_start

            if this_end < last_end:
                # this is contained

                if this_score > last_score:
                    # print initial fragment of last
                    init_a = copy.copy(last_a)
                    init_a[6] = str(this_start-1)
                    print gff_line(init_a)

                    # print this
                    print gff_line(this_a)

                    # save final fragment of last
                    last_a[5] =  str(this_end+1)

                else:
                    # skip this
                    pass

            else:
                # this overlaps to the right
                if this_score > last_score:
                    # print fragment of last up to this
                    last_a[6] = str(this_start-1)
                    print gff_line(last_a)

                    # save full this
                    last_a = this_a
                    
                else:
                    # print full fragment of last
                    print gff_line(last_a)

                    # save fragment of this
                    this_a[5] = str(last_end+1)
                    last_a = this_a

        line = rm_in.readline()

    print gff_line(last_a)


################################################################################
# gff_line
################################################################################
def gff_line(a):
    strand = a[8]
    if strand == 'C':
        strand = '-'

    cols = (a[4], 'RepeatMasker', 'repeat', a[5], a[6], '.', strand, '.', gff.kv_gtf({'repeat':a[9], 'family':a[10]}))
    return '\t'.join(cols)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
