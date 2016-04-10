#!/usr/bin/env python
from optparse import OptionParser

################################################################################
# possum2bed.py
#
#
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <possum_file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide possum output file')
    else:
        possum_file = args[0]

    for line in open(possum_file):
        a = line.split('\t')

        tf_id = a[0]
        start = int(a[5])+1
        end = start+int(a[6])-1
        fnrc = a[7]
        score = a[9]
        seq_id = a[16][:a[16].find('.')]

        if fnrc == 'fn':
            strand = '+'
        else:
            strand = '-'

        out_a = [seq_id, str(start), str(end), tf_id, score, strand]
        print '\t'.join(out_a)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
