#!/usr/bin/env python
from optparse import OptionParser
import pdb, sys
import gff

################################################################################
# tgff_cgff.py
#
# Map transcript features to chromosomes when the transcripts can be spliced.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <transcript .gff>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='cgff_file', default='/Users/dk/research/common/data/lncrna/lnc_catalog.gtf', help='Gtf file mapping transcripts to chromosomes [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide gff file mapping features to transcripts')
    else:
        tgff_file = args[0]

    # get transcript information
    transcripts = {}
    for line in open(options.cgff_file):
        a = line.split('\t')
        if a[2] == 'exon':
            trans_id = gff.gtf_kv(a[8])['transcript_id']
            if not trans_id in transcripts:
                transcripts[trans_id] = Transcript(trans_id,a[0],a[6])
            transcripts[trans_id].add_exon(int(a[3]), int(a[4]))

    # process transcript features
    for line in open(tgff_file):
        feat = Feature(line)
        map_feature(transcripts[feat.trans_id], feat)


################################################################################
# flip_feature
################################################################################
def flip_feature(transcript, feat):
    trans_length = 0
    for (estart,eend) in transcript.exons:
        trans_length += eend-estart+1

    feat_length = feat.end - feat.start

    feat.start = trans_length - feat.end + 1
    feat.end = feat.start + feat_length
    if feat.strand == '+':
        feat.strand = '-'
    else:
        feat.strand = '+'


################################################################################
# map_feature
################################################################################
def map_feature(transcript, feat):    
    if transcript.strand == '-':
        # flip feature and treat as '+'
        flip_feature(transcript, feat)
    
    events = [Event(feat.start, -1, 'fstart'), Event(feat.end, -1, 'fend')]
    tstart = transcript.exons[0]
    tbases = 0
    for (estart,eend) in transcript.exons:
        events.append(Event(tbases+1, estart, 'estart'))
        events.append(Event(tbases+1+eend-estart, eend, 'eend'))
        tbases += eend-estart+1
    events[2].type = 'tstart'
    events[-1].type = 'tend'

    '''
    print >> sys.stderr, transcript.id
    for e in events:
        print >> sys.stderr, e
    print >> sys.stderr, ''
    '''

    events.sort()

    in_feature = False
    for i in range(len(events)-1):
        if events[i].type == 'fstart':
            # position feature start
            in_feature = True
            events[i].cpos = events[i-1].cpos + events[i].tpos - events[i-1].tpos

        elif in_feature:            
            if events[i].type == 'fend':
                # finish feature
                events[i].cpos = events[i-1].cpos + events[i].tpos - events[i-1].tpos
                cols = [transcript.chrom, feat.source, feat.name, str(events[i-1].cpos), str(events[i].cpos), '.', feat.strand, '.', feat.group+' '+transcript.id]
                print '\t'.join(cols)
                break

            elif events[i].type == 'eend':
                # feature is spliced
                cols = [transcript.chrom, feat.source, feat.name, str(events[i-1].cpos), str(events[i].cpos), '.', feat.strand, '.', feat.group+' '+transcript.id]
                print '\t'.join(cols)

            elif events[i].type == 'estart':
                pass

            else:
                print >> sys.stderr, 'Unimagined event case'
                exit(1)

class Feature:
    def __init__(self, gff_line):
        a = gff_line.split('\t')
        self.trans_id = a[0]
        self.source = a[1]
        self.name = a[2]
        self.start = int(a[3])
        self.end = int(a[4])
        self.strand = a[6]
        self.group = a[8].rstrip()

class Transcript:
    def __init__(self, tid, chrom, strand):
        self.id = tid
        self.chrom = chrom
        self.strand = strand
        self.exons = []
    def add_exon(self, start, end):
        self.exons.append((start,end))

class Event:
    def __init__(self, tpos, cpos, etype):
        self.type = etype
        self.cpos = cpos
        self.tpos = tpos

    def __str__(self):
        return '%9d %9d %8s' % (self.tpos, self.cpos, self.type)

    def __cmp__(self, other):
        if self.tpos == other.tpos:
            if self.type[1:] == 'start' and other.type[1:] == 'end':
                return -1
            elif self.type[1:] == 'end' and other.type[1:] == 'start':
                return 1

            elif self.type == 'tstart' and other.type == 'fstart':
                return -1
            elif self.type == 'fstart' and other.type == 'tstart':
                return 1

            elif self.type == 'estart' and other.type == 'fstart':
                return -1
            elif self.type == 'fstart' and other.type == 'estart':
                return 1

            elif self.type == 'fend' and other.type == 'tend':
                return -1
            elif self.type == 'tend' and other.type == 'fend':
                return 1

            elif self.type == 'fend' and other.type == 'eend':
                return -1
            elif self.type == 'eend' and other.type == 'fend':
                return 1

            else:
                print >> sys.stderr, 'Unimagined sort case'
                print >> sys.stderr, self
                print >> sys.stderr, other
                exit(1)

        return self.tpos - other.tpos

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
