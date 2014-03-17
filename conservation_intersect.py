#!/usr/bin/env python
from optparse import OptionParser
from bx.intervals.intersection import Interval, IntervalTree
import gzip, glob, os, sys, subprocess

################################################################################
# conservation_intersect.py
#
# Intersect a list of segments (e.g. lincRNAs) in gff format with the multiZ
# blocks and print out the phastCons/PhyloP scores.
#
# Assumes that the gff entries are disjoint which can be accomplished using
# mergeBed.
#
# mergeBed has a little quirk where a 1 bp gff entry will be changed into
# a 2 bp entry, which causes very slight differences between using the '-l'
# option and not.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gff file>'
    parser = OptionParser(usage)
    parser.add_option('-b', dest='background_gff', help='GFF file describing valid background sequences [Default: %default] NOT IMPLEMENTED')
    parser.add_option('-l', dest='lncrna', action='store_true', default=False, help='Use the lncRNA specific phastcons file to speed things up [Default: %default]')
    parser.add_option('-c', dest='conservation_type', default='phylop', help='Conservation type to use [phastcons|phylop] [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide gff file to intersect')
    gff_file = args[0]

    cons_dir = '%s/research/common/data/%s' % (os.environ['HOME'],options.conservation_type)
    if not os.path.isdir(cons_dir):
        parser.error('Must specify conservation type as "phylop" or "phastcons"')    

    # build interval trees
    print >> sys.stderr, 'Building interval trees ...',
    chr_features = {}
    p = subprocess.Popen('sortBed -i %s | mergeBed -i -' % gff_file, shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        a = line.split('\t')
        chr_features.setdefault(a[0], IntervalTree()).insert_interval( Interval(int(a[1])+1,int(a[2])) )
    p.communicate()
    print >> sys.stderr, 'Done'

    # build background interval trees
    if options.background_gff:
        chr_features_bg = sample_background_intervals(gff_file, options.background_gff)
        
    # process overlapping chromosome blocks
    if options.lncrna:
        lnc_files = glob.glob('%s/lnc_catalog.*wigFix.gz' % cons_dir)
        if len(lnc_files) != 1:
            print >> sys.stderr, 'Ambiguous lnc catalog file'
            exit(1)
        else:
            process_file(chr_features, lnc_files[0])

    else:
        for pc_file in glob.glob('%s/chr*' % cons_dir):
            process_file(chr_features, pc_file)


################################################################################
# intersect_scores
#
# Print out block scores overlapping features.
################################################################################
def intersect_scores(features, block_start, block_scores):
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

        print '\n'.join([str(s) for s in block_scores[start:end]])


################################################################################
# process_file
#
# Process overlapping chromosome blocks in the given file.
################################################################################
def process_file(chr_features, pc_file):
    if pc_file[-2:] == 'gz':
        pc_f = gzip.open(pc_file)
    elif os.path.isfile(pc_file):
        pc_f = open(pc_file)
    elif os.path.isfile(pc_file+'.gz'):
        pc_f = gzip.open(pc_file+'.gz')

    chrom = os.path.split(pc_file)[1].split('.')[0]
    print >> sys.stderr, 'Processing %s ...' % chrom,

    block_start = 0
    block_scores = []

    line = pc_f.readline()
    while line:
        if line.startswith('fixedStep'):
            if block_scores:
                intersect_scores(chr_features.get(chrom, IntervalTree()), block_start, block_scores)

            a = line.split()
            chrom = a[1][6:]
            block_start = int(a[2][6:])
            block_scores = []
        else:
            block_scores.append(float(line.rstrip()))

        line = pc_f.readline()

    intersect_scores(chr_features.get(chrom, IntervalTree()), block_start, block_scores)

    pc_f.close()
    print >> sys.stderr, 'Done'


################################################################################
# sample_background_intervals
#
# Simulate many datasets of the same size and length distribution from the
# background sequence specificed.
#
# Input
#  feature_gff:    GFF file of features.
#  background_gff: GFF file of background sequences from which to simulate
#                   features.
#
# Output
#  ???
################################################################################
def sample_background_intervals(feature_gff, background_gff):
    # determine length distributions
    feature_lengths = {}
    for line in open(feature_gff):
        a = line.split('\t')
        flen = int(a[4]) - int(a[3]) + 1
        feature_lengths[flen] = feature_lengths.get(flen,0) + 1

    # merge background gff
    background_intervals = []
    p = subprocess.Popen('sortBed -i %s | mergeBed -i -' % background_gff, shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        a = line.split('\t')
        chrom = a[0]
        start = int(a[1])
        end = int(a[2])

        # convert to gff
        background_intervals.append((chrom,start+1,end))

    # sample
    chr_features_bg = {}


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
