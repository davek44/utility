#!/usr/bin/env python
from optparse import OptionParser
import gzip, sys

############################################################
# cutFasta
#
# Extract a sequence from a fasta file, using 1-based
# indexing
############################################################


############################################################
# main
############################################################
def main():
    parser = OptionParser()
    parser.add_option('-x', dest='start', type='int', help='Cut start')
    parser.add_option('-y', dest='end', type='int', help='Cut end')
    parser.add_option('-s', dest='header', help='Sequence header')
    parser.add_option('-c', dest='contain', action='store_true', default=False, help='Grab all sequences that contain the header pattern')
    (options,args) = parser.parse_args()
    
    if len(args) != 1:
        parser.error('Please provide a single fasta file')
        
    print cf(options.start, options.end, options.header, options.contain, args[0])

############################################################
# cf
#
# Pull out the sequence from start to end in the entry
# header in the file fasta_file
############################################################
def cf(start, end, header, contain, fasta_file):
    # collect sequence up to end
    seq = ''
    get_seq = False

    if fasta_file[-3:] == '.gz':
        ff = gzip.open(fasta_file)
    else:
        ff = open(fasta_file)
    line = ff.readline()
    while line:
        if line[0] == '>':
            if get_seq:
                # already found, stop
                break
            else:
                # check header
                h = line[1:].rstrip()
                if not header:
                    get_seq = True
                    header = h
                elif h == header or (contain and h.find(header) != -1):
                    get_seq = True

        elif get_seq:
            seq += line.rstrip()

        # if past end, stop
        if end and len(seq) > end:
            break

        line = ff.readline()

    # print seq
    if start and end:
        return '>%s_(%d-%d)\n%s' % (header,start,end,seq[start-1:end])
    elif start:
        return '>%s_(%d-)\n%s' % (header,start,seq[start-1:])
    elif end:
        return '>%s_(-%d)\n%s' % (header,start,seq[:end])
    else:
        return '>%s\n%s' % (header,seq)

############################################################
# __main__
############################################################
if __name__ == '__main__':
    main()
