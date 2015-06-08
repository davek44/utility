#!/usr/bin/env python
from optparse import OptionParser
from bx.intervals.intersection import Interval, IntervalTree
import gzip, glob, os, random, shutil, sys, subprocess, tempfile

################################################################################
# conservation_intersect_motif.py
#
# Intersect a list of motifs in GFF format with the multiZ blocks and print out
# the phastCons/PhyloP scores for each motif position.
#
# Assumes that the gff entries are disjoint.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gff file>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='conservation_type', default='phylop', help='Conservation type to use [phastcons|phylop] [Default: %default]')
    parser.add_option('-m', dest='max_instances', default=5000, type='int', help='Maximum number of instances to consider per TE-motif [Default: %default]')
    parser.add_option('-p', dest='position', default=False, action='store_true', help='Print the positional index of the conservation value within the feature. (May allow duplicated values). [Default: %default]')
    parser.add_option('-o', dest='out_dir', help='Output directory')
    parser.add_option('-r', dest='range', type='int', default=0, help='Range to expand each motif occurrence to [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide gff file to intersect')
    gff_file = args[0]

    cons_dir = '%s/research/common/data/%s' % (os.environ['HOME'],options.conservation_type)
    if not os.path.isdir(cons_dir):
        parser.error('Must specify conservation type as "phylop" or "phastcons"')

    ############################################################
    # determine sample prob
    ############################################################
    # count instances per motif
    motif_counts = {}
    for line in open(gff_file):
        a = line.split('\t')
        motif_id = a[8]
        motif_counts[motif_id] = motif_counts.get(motif_id,0) + 1

    # compute prob
    motif_sample = {}
    for motif_id in motif_counts:
        motif_sample[motif_id] = min(1.0, options.max_instances / float(motif_counts[motif_id]))

    ############################################################
    # extend GFF entries to range (and sample)
    ############################################################
    if options.range > 0:
        range_gff_fd, range_gff_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
        range_gff_out = open(range_gff_file, 'w')
        for line in open(gff_file):
            a = line.split('\t')

            pstart = int(a[3])
            pend = int(a[4])
            motif_id = a[8]

            if random.random() < motif_sample[motif_id]:
                peak_mid = pstart + (pend-pstart)/2
                a[3] = str(peak_mid - options.range/2)
                a[4] = str(peak_mid + options.range/2)

                print >> range_gff_out, '\t'.join(a),
        range_gff_out.close()
    else:
        range_gff_file = gff_file

    ############################################################
    # determine where to get conservation values
    ############################################################
    # initializing data structures
    motif_chr_features = {}
    for line in open(range_gff_file):
        a = line.split('\t')
        motif_id = a[8].rstrip()
        motif_chr_features[motif_id] = {}

    # build interval trees
    print >> sys.stderr, 'Building interval trees ...'
    if options.position:
        for line in open(range_gff_file):
            a = line.split('\t')
            motif_id = a[8].rstrip()
            motif_chr_features[motif_id].setdefault(a[0], IntervalTree()).insert_interval( Interval(int(a[3]),int(a[4])) )
    else:
        for motif_id in motif_chr_features:
            p = subprocess.Popen('awk \'$9 == "%s"\' %s | sortBed -i - | mergeBed -i -' % (motif_id,range_gff_file), shell=True, stdout=subprocess.PIPE)
            for line in p.stdout:
                a = line.split('\t')
                motif_chr_features[motif_id].setdefault(a[0], IntervalTree()).insert_interval( Interval(int(a[1])+1,int(a[2])) )
            p.communicate()
    print >> sys.stderr, 'Done'

    ############################################################
    # get conservation values
    ############################################################
    # open output files
    if os.path.isdir(options.out_dir):
        shutil.rmtree(options.out_dir)
    os.mkdir(options.out_dir)
    
    motif_outs = {}
    for motif_id in motif_chr_features:
        motif_outs[motif_id] = open('%s/%s.txt' % (options.out_dir,motif_id), 'w')
        
    # process overlapping chromosome blocks
    for pc_file in glob.glob('%s/chr*' % cons_dir):
        process_file(motif_chr_features, pc_file, motif_outs, options.position)

    # close output files
    for motif_id in motif_outs:
        motif_outs[motif_id].close()

    # clean
    if options.range > 0:
        os.close(range_gff_fd)
        os.remove(range_gff_file)


################################################################################
# intersect_scores
#
# Print out block scores overlapping features.
################################################################################
def intersect_scores(features, block_start, block_scores, out_open, position):
    block_end = block_start+len(block_scores)-1
    for overlap_interval in features.find(block_start, block_start+len(block_scores)):
        # block internal to interval
        if overlap_interval.start <= block_start <= block_end <= overlap_interval.end:
            start = 0
            end = len(block_scores)
            feature_start = block_start - overlap_interval.start

        # interval internal to block
        elif block_start <= overlap_interval.start <= overlap_interval.end <= block_end:
            start = overlap_interval.start - block_start
            end = start + overlap_interval.end - overlap_interval.start + 1
            feature_start = 0

        # left block overlap interval
        elif block_start < overlap_interval.start:
            start = overlap_interval.start - block_start
            end = start + block_end - overlap_interval.start + 1
            feature_start = 0
            
        # right block overlap interval
        else:
            start = 0
            end = overlap_interval.end - block_start + 1
            feature_start = block_start - overlap_interval.start

        if not position:
            print >> out_open, '\n'.join([str(s) for s in block_scores[start:end]])
        else:
            olap_len = end - start
            for i in range(olap_len):
                print >> out_open, block_scores[start+i], (feature_start+i)


################################################################################
# process_file
#
# Process overlapping chromosome blocks in the given file.
################################################################################
def process_file(motif_chr_features, pc_file, motif_outs, position):
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
                for motif_id in motif_chr_features:
                    intersect_scores(motif_chr_features[motif_id].get(chrom, IntervalTree()), block_start, block_scores, motif_outs[motif_id], position)

            a = line.split()
            chrom = a[1][6:]
            block_start = int(a[2][6:])
            block_scores = []
        else:
            block_scores.append(float(line.rstrip()))

        line = pc_f.readline()

    for motif_id in motif_chr_features:
        intersect_scores(motif_chr_features[motif_id].get(chrom, IntervalTree()), block_start, block_scores, motif_outs[motif_id], position)

    pc_f.close()
    print >> sys.stderr, 'Done'


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
