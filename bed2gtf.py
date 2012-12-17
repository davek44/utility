#!/usr/bin/env python
from optparse import OptionParser

################################################################################
# bed2gtf.py
#
# Convert a bed file to a gtf file.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide gtf file')
    else:
        bed_file = args[0]

    for line in open(bed_file):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        tid = a[3]

        gene_start = int(a[1])
        gene_end = int(a[2])

        block_sizes = [int(x) for x in a[10].split(',') if x]
        block_starts = [int(x) for x in a[11].split(',') if x]

        exon_num = 1
        for i in range(len(block_starts)):
            exon_start = gene_start+1+block_starts[i]
            exon_end = gene_start+1+block_starts[i]+block_sizes[i]-1

            cols = [a[0], 'BED', 'exon', str(exon_start), str(exon_end), '.', a[5], '.', 'transcript_id "%s"; exon_number "%d"' % (tid,exon_num)]
            print '\t'.join(cols)
            exon_num += 1
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
