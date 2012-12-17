#!/usr/bin/env python
from optparse import OptionParser
from bx.intervals.intersection import Interval, IntervalTree
import gzip, pdb, sys
import gff

################################################################################
# chain_map.py
#
# Map annotations in gff format to a new genome via a chain file.
#
# Assuming that features with the same group attribute are all on the same
# strand and can be treated as one big feature.
#
# The chain files are messy and contain overlapping alignment blocks, so
# I'm separating each chain. To further filter the output, provide a net
# file, and I'll only map to chains in the net.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <chain file> <gff file>'
    parser = OptionParser(usage)
    parser.add_option('-k', dest='gtf_key', help='Group based on the given gtf key [Default: %default]')
    parser.add_option('-m', dest='merge_t', type='int', default=40, help='Minimum distance between alignment blocks to merge into a single exon [Default: %default]')
    parser.add_option('-n', dest='net_file', help='Net file to filter the chains considered')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error(usage)
    else:
        chain_file = args[0]
        gff_file = args[1]

    # build interval tree's of gff annotations
    # map intervals to annotations
    chr_features, interval_map = gff_intervals(gff_file, options.gtf_key)

    # filter by net chains
    if options.net_file:
        net_chains = get_net_chains(options.net_file)
    else:
        net_chains = None

    # process chain file
    if chain_file[-3:] == '.gz':
        chain_in = gzip.open(chain_file)
    else:
        chain_in = open(chain_file)

    # header
    line = chain_in.readline()
    while line[:2] == '##':
        line = chain_in.readline()

    # map features
    mapped_features = {}
    chain_lines = []
    while line:
        if line.startswith('chain'):
            if chain_lines and (not net_chains or chain_id in net_chains):
                process_chain(chr_features, interval_map, mapped_features, chain_lines)
            chain_lines = [line]
            chain_id = line.split()[-1]
        elif line.rstrip(): # skip blank lines
            chain_lines.append(line)

        line = chain_in.readline()
    chain_in.close()

    if chain_lines and (not net_chains or chain_id in net_chains):
        process_chain(chr_features, interval_map, mapped_features, chain_lines)

    # merge mapped features as gff lines
    for feature_id in mapped_features:
        for chain_id in mapped_features[feature_id]:
            # sort features
            mapped_features[feature_id][chain_id].sort(feature_cmp)

            # make first line
            (qchrom,qstart,qend,qstrand) = mapped_features[feature_id][chain_id][0]
            gff_cols = [[qchrom, 'MyTransMap', 'feature', qstart, qend, '.', qstrand, '.', 'chain_id %s; %s'%(chain_id,feature_id)]]

            # make the rest of the lines
            for (qchrom,qstart,qend,qstrand) in mapped_features[feature_id][chain_id][1:]:

                if gff_cols[-1][0] == qchrom and gff_cols[-1][6] == qstrand and gff_cols[-1][4] + options.merge_t >= qstart:
                    # large overlaps unexpected
                    ovl_start = max(gff_cols[-1][3], qstart)
                    ovl_end = min(gff_cols[-1][4], qend)
                    if (ovl_end-ovl_start+1) > 50:
                        print >> sys.stderr, 'Large overlap between features - %s' % feature_id

                    # merge w/ prior
                    gff_cols[-1][4] = qend
                else:
                    # add new
                    gff_cols.append([qchrom, 'MyTransMap', 'feature', qstart, qend, '.', qstrand, '.', 'chain_id %s; %s'%(chain_id,feature_id)])

            # print lines
            for cols in gff_cols:
                print '\t'.join([str(c) for c in cols])


################################################################################
# get_net_chains
#
# Save all of the chains used in nets to a set.
################################################################################
def get_net_chains(net_file):
    net_chains = set()

    if net_file[-3:] == '.gz':
        net_in = gzip.open(net_file)
    else:
        net_in = open(net_file)

    line = net_in.readline()
    while line:
        a = line.split()
        if a[0] == 'fill':
            for i in range(len(a)):
                if a[i] == 'id':
                    net_chains.add(a[i+1])
                    break
        line = net_in.readline()
                
    return net_chains


################################################################################
# gff_intervals
#
# Build interval tree's of gff annotations and map intervals to annotations.
################################################################################
def gff_intervals(gff_file, gtf_key):
    chr_features = {}
    interval_map = {}

    for line in open(gff_file):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        chrom = a[0]
        start = int(a[3])
        end = int(a[4])
        strand = a[6]
        if gtf_key:
            feature_id = gff.gtf_kv(a[8]).get(gtf_key,a[8])
        else:
            feature_id = a[8]

        chr_features.setdefault(chrom,IntervalTree()).insert_interval(Interval(start,end))
        interval_map.setdefault(chrom,{}).setdefault((start,end),[]).append((feature_id,strand))

    return chr_features, interval_map


################################################################################
# feature_cmp
#
# Sort mapped features before merging.
################################################################################
def feature_cmp(x, y):
    # chrom
    if x[0] < y[0]:
        return -1
    elif x[0] > y[0]:
        return 1
    else:
        # strand
        if x[3] < y[3]:
            return -1
        elif x[3] > y[3]:
            return 1
        else:
            # start
            if x[1] < y[1]:
                return -1
            elif x[1] > y[1]:
                return 1
            else:
                # end
                if x[2] < y[2]:
                    return -1
                else:
                    return 1

    
################################################################################
# process_chain
#
# Map the annotation features described in chr_features and interval_map
# from the reference to query genomes described by the alignment chain in
# chain_lines and save the mappings to mapped_features.
################################################################################
def process_chain(chr_features, interval_map, mapped_features, chain_lines):
    a = chain_lines[0].split()
    rchrom = a[2]
    rstrand = a[4]
    rstart = int(a[5])
    rend = int(a[6])
    qchrom = a[7]
    qstrand = a[9]
    if qstrand == '+':
        qstart = int(a[10])
        qend = int(a[11])
    else:
        qstart = int(a[8]) - int(a[11])
        qend = int(a[8]) - int(a[10])
    chain_id = a[-1]
    map_dir = (rstrand == qstrand)

    # to follow along w/ alignments in the chain
    align_rstart = rstart
    if map_dir:
        align_qstart = qstart
    else:
        align_qend = qend

    if rchrom in chr_features:
        for line in chain_lines[1:]:            
            # get alignment info
            a = line.split()
            align_size = int(a[0])
            if len(a) > 1:
                rgap = int(a[1])
                qgap = int(a[2])
            else:
                rgap = 0
                qgap = 0

            align_rend = align_rstart + align_size
            if map_dir:
                align_qend = align_qstart + align_size
            else:
                align_qstart = align_qend - align_size

            # search for overlapping feature intervals
            for feature_interval in chr_features[rchrom].find(align_rstart, align_rend):
                # can only map the overlapping portion
                map_rstart = max(align_rstart, feature_interval.start)
                map_rend = min(align_rend, feature_interval.end)

                # map                
                if map_dir:
                    map_qstart = align_qstart + (map_rstart - align_rstart)
                    map_qend = align_qstart + (map_rend - align_rstart)
                else:
                    map_qstart = align_qstart + (align_rend - map_rend)
                    map_qend = align_qstart + (align_rend - map_rstart)
                
                #log_cols = [chain_id, align_rstart, align_rend, align_qstart, align_qend, feature_interval.start, feature_interval.end, map_rstart, map_rend, map_qstart, map_qend]
                #print >> sys.stderr, '\t'.join([str(c) for c in log_cols])

                # apply to all features w/ this interval
                for (feature_id,map_rstrand) in interval_map[rchrom][(feature_interval.start,feature_interval.end)]:
                    map_qstrand = map_rstrand
                    if not map_dir:
                        # flip
                        if map_qstrand == '+':
                            map_qstrand = '-'
                        else:
                            map_qstrand = '+'

                    # add 1 to start to convert 0-based to 1-based
                    mapped_features.setdefault(feature_id,{}).setdefault(chain_id,[]).append((qchrom,map_qstart+1,map_qend,map_qstrand))

            # get past the gap
            align_rstart = align_rend + rgap
            if map_dir:
                align_qstart = align_qend + qgap
            else:
                align_qend = align_qstart - qgap


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
