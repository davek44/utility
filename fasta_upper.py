#!/usr/bin/env python
from optparse import OptionParser

'''
Name

Description...
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_fasta_file> <out_fasta_file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    in_fasta_file = args[0]
    out_fasta_file = args[1]
    out_fasta_open = open(out_fasta_file, 'w')

    for line in open(in_fasta_file):
        if line[0] == '>':
            print(line, end='', file=out_fasta_open)
        else:
            print(line.upper(), end='', file=out_fasta_open)
            
    out_fasta_open.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
