#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess, random, sys

################################################################################
# blastme.py
#
# My wrapper for Blast so that I can ease into it from Mummer.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <subject> <query>'
    parser = OptionParser(usage)
    parser.add_option('-b', dest='blast', default='blastn', help='Blast program [Default: %default]')
    parser.add_option('-c', dest='cover', default=0, type='float', help='Coverage% of query sequence required to output alignment [Default: %default]')
    parser.add_option('-e', '--evalue', dest='evalue', type='float', default=1e-5, help='E-value threshold [Default: %default]')
    parser.add_option('-i', '--perc_identity', dest='perc_identity', type='float', default=0, help='Percent identity threshold [Default: %default]')
    parser.add_option('-o', dest='opt_str', help='Options string to pass to blast')
    parser.add_option('-p', '--num_threads', dest='num_threads', type='int', default=1, help='Index word size [Default: %default]. GIVES STOCHASTIC RESULTS!')
    parser.add_option('-s', '--self', dest='self_align', action='store_true', default=False, help='Self alignment; filter output accordingly [Default: %default]')
    parser.add_option('-w', '--word_size', dest='word_size', type='int', default=28, help='Index word size [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide subject and query files')
    else:
        subject_file = args[0]
        query_file = args[1]

    if options.num_threads != 1:
        print >> sys.stderr, 'Warning: multi-threaded programs seem to reduce sensitivity and add stochasticity'

    # get subject names, lengths
    subject_names = []
    subject_lengths = {}
    for line in open(subject_file):
        if line[0] == '>':
            header = line[1:].split()[0]
            subject_names.append(header)
            if header in subject_lengths:
                print >> sys.stderr, 'Double header - %s' % header
            subject_lengths[header] = 0
        else:
            subject_lengths[header] += len(line.rstrip())

    # get querylengths
    query_lengths = {}
    for line in open(query_file):
        if line[0] == '>':
            header = line[1:].split()[0]
            if header in query_lengths:
                print >> sys.stderr, 'Double header - %s' % header
            query_lengths[header] = 0
        else:
            query_lengths[header] += len(line.rstrip())
    
    # BLAST
    r = random.randint(0,999999)
    cmd = '%s -subject %s -query %s -evalue %f -perc_identity %f -num_threads %d -word_size %d -outfmt 6 > %d.tmp' % (options.blast, subject_file, query_file, options.evalue, options.perc_identity, options.num_threads, options.word_size, r)
    print >> sys.stderr, cmd
    p = subprocess.Popen(cmd, shell=True)
    os.waitpid(p.pid, 0)

    # parse output
    for line in open('%d.tmp' % r):
        a = line.split()
        
        stag = subject_names[int(a[1][8:])-1]
        qtag = a[0]

        if options.self_align:
            if stag >= qtag:
                continue

        sstart = int(a[8])
        send = int(a[9])
        qstart = int(a[6])
        qend = int(a[7])

        idy = float(a[2])/100.0
        evalue = a[10]

        if sstart > send:
            sstart, send = send, sstart
            qstart, qend = qend, qstart

        if float(abs(qend-qstart)+1) / query_lengths[qtag] > options.cover:
            print '%9d %9d  | %9d %9d  | %8d %8d | %9d %9d  |  %.4f %7s  | %9s %9s' % (sstart, send, qstart, qend, (send-sstart+1), (abs(qend-qstart)+1), subject_lengths[stag], query_lengths[qtag], idy, evalue, stag, qtag)

    os.remove('%d.tmp' % r)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
