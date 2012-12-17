#!/usr/bin/env python
from optparse import OptionParser
from bx.intervals.intersection import Interval, IntervalTree
import gzip, glob, os, sys
import gff, stats

################################################################################
# lnc_phylop.py
#
# Intersect a set of lncRNAs in gtf format with the multiZ blocks and compute
# stats about the the PhyloP scores.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gff file>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='cons_dir', default='%s/research/common/data/phylop' % os.environ['HOME'], help='Conservation directory [Default: %default]')
    parser.add_option('-l', dest='lncrna', action='store_true', default=False, help='Use the lncRNA specific file to speed things up [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide gff file to intersect')
    else:
        gff_file = args[0]

    t2g = gff.t2g(gff_file)

    # build interval trees
    lnc_lengths = {}
    chr_features = {}
    interval2lnc = {}
    lnc_cons = {}
    for line in open(gff_file):
        a = line.split('\t')

        chrom = a[0]
        start = int(a[3])
        end = int(a[4])
        tid = gff.gtf_kv(a[8])['transcript_id']
        align = (chrom,start,end)

        lnc_cons[tid] = []
        lnc_lengths[tid] = lnc_lengths.get(tid,0) + (end-start+1)
        if interval2lnc.has_key(align):
            interval2lnc[align].add(tid)
        else:
            interval2lnc[align] = set([tid])
            chr_features.setdefault(chrom, IntervalTree()).insert_interval(Interval(start,end))

    # process overlapping chromosome blocks
    if options.lncrna:
        lnc_wig = glob.glob('%s/lnc_catalog.*wigFix*' % options.cons_dir)[0]
        process_file(chr_features, interval2lnc, lnc_cons, lnc_wig)

    else:
        for cons_file in glob.glob('%s/chr*' % options.cons_dir):
            process_file(chr_features, interval2lnc, lnc_cons, cons_file)

    # print table
    for tid in lnc_lengths:
        cons_len = len(lnc_cons[tid])
        cons_cov = float(cons_len) / lnc_lengths[tid]
        if cons_len == 0:
            cons_mean = 0.0
            cons_median = 0.0
            cons_pos = 0.0
            cons_neg = 0.0
        else:
            cons_mean = stats.mean(lnc_cons[tid])
            cons_median = stats.median(lnc_cons[tid])
            cons_pos = len([c for c in lnc_cons[tid] if c > 1]) / float(cons_len)
            cons_neg = len([c for c in lnc_cons[tid] if c < 1]) / float(cons_len)

        cols = (tid, t2g[tid], lnc_lengths[tid], cons_cov, cons_mean, cons_median, cons_neg, cons_pos)
        print '%-15s %-15s %7d %9.4f %9.4f %9.4f %9.4f %9.4f' % cols


################################################################################
# intersect_scores
#
# Print out block scores overlapping features.
################################################################################
def intersect_scores(chr_features, interval2lnc, lnc_cons, chrom, block_start, block_scores):
    features = chr_features.get(chrom, IntervalTree())
    block_end = block_start+len(block_scores)-1
    for overlap_interval in features.find(block_start, block_start+len(block_scores)):
        # block internal to invterval
        if overlap_interval.start <= block_start <= block_end <= overlap_interval.end:
            start = 0
            end = len(block_scores)

        # interval internal to block
        elif block_start <= overlap_interval.start <= overlap_interval.end <= block_end:
            start = overlap_interval.start - block_start
            end = start + overlap_interval.end - overlap_interval.start + 1

        # left block overlap interval
        elif block_start < overlap_interval.start:
            start = overlap_interval.start - block_start
            end = start + block_end - overlap_interval.start + 1
            
        # right block overlap interval
        else:
            start = 0
            end = overlap_interval.end - block_start + 1

        #start = overlap_interval.start - block_start
        #end = start + overlap_interval.end - overlap_interval.start

        for tid in interval2lnc[(chrom,overlap_interval.start,overlap_interval.end)]:
            lnc_cons[tid] += block_scores[start:end]


################################################################################
# process_file
#
# Process overlapping chromosome blocks in the given file.
################################################################################
def process_file(chr_features, interval2lnc, lnc_cons, cons_file):        
    if cons_file[-2:] == 'gz':
        cons_f = gzip.open(cons_file)
    else:
        cons_f = open(cons_file)

    chrom = os.path.split(cons_file)[1].split('.')[0]
    print >> sys.stderr, 'Processing %s ...' % chrom,

    block_start = 0
    block_scores = []

    line = cons_f.readline()
    while line:
        if line.startswith('fixedStep'):
            if block_scores:
                intersect_scores(chr_features, interval2lnc, lnc_cons, chrom, block_start, block_scores)

            a = line.split()
            chrom = a[1][6:]            
            block_start = int(a[2][6:])
            block_scores = []
        else:
            block_scores.append(float(line.rstrip()))

        line = cons_f.readline()

    intersect_scores(chr_features, interval2lnc, lnc_cons, chrom, block_start, block_scores)

    cons_f.close()
    print >> sys.stderr, 'Done'



################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
