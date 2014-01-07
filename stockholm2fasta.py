#!/usr/bin/env python
from optparse import OptionParser

################################################################################
# stockholm2fasta.py
#
# Convert Stockholm MSA format from HMMer to FASTA for viewing.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <stockholm>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='consensus_only', default=False, action='store_true', help='Print consensus columns only [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide input Stockholm format file')
    else:
        stockholm_file = args[0]

    seqs = {}
    consensus = ''
    for line in open(stockholm_file):
        if line.rstrip() not in ['','//'] and line[0] != '#':
            header, msa = line.split()
            seqs[header] = seqs.get(header,'') + msa
        elif line.startswith('#=GC RF'):
            consensus += line.split()[-1]

    if options.consensus_only:
        for header in seqs:
            hseq = seqs[header]
            seqs[header] = ''.join([hseq[i] for i in range(len(hseq)) if consensus[i] == 'x'])

    for header in seqs:
        print '>%s\n%s' % (header,seqs[header])


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
