#!/usr/bin/env python
from optparse import OptionParser
import re

################################################################################
# trf_mask.py
#
# Mask the tandem repeats found by Tandem Repeat Finder.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <seq file> <trf file 1> ... <trf file N>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) < 2:
        parser.error('Please provide sequence file and TRF output file')
    else:
        seq_file = args[0]
        trf_files = args[1:]

    repeats = {}
    for trf_file in trf_files:
        get_repeats(trf_file, repeats)

    header = ''
    for line in open(seq_file):
        if line[0] == '>':
            if header:
                mseq = mask_seq(seq, repeats.get(header,[]))
                print '>%s\n%s' % (header,mseq)

            header = line[1:].rstrip()
            seq = ''

        else:
            seq += line.rstrip()

    if header:
        mseq = mask_seq(seq, repeats.get(header,[]))
        print '>%s\n%s' % (header,mseq)


################################################################################
# get_repeats
#
# Save the repeats in a dict keyed by the header
################################################################################
def get_repeats(trf_file, repeats):
    indices_re = re.compile('Indices: (\d+)\-\-(\d+)\s*Score')
    for line in open(trf_file):
        if line.startswith('Sequence:'):
            header = line[10:].rstrip()
        else:
            m = indices_re.search(line)
            if m:
                (start,end) = m.group(1,2)
                repeats.setdefault(header,[]).append((int(start)-1,int(end)))


################################################################################
# mask_seq
#
# Mask the sequence using the list of repeats
################################################################################
def mask_seq(seq, seq_repeats):
    mseq = list(seq)
    for rep in seq_repeats:
        for i in range(rep[0],rep[1]):
            mseq[i] = 'N'
    return ''.join(mseq)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
