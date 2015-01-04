#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import binom, norm
from te_bam_enrich import te_target_size, count_hg19, count_bed
import gzip, os, subprocess, sys, tempfile
import fdr, gff, stats

################################################################################
# te_gff_enrich.py
#
# Compute statistics about the overlap occurrences of GFF entries on TEs.
#
# Note:
#  -The statistics here assume consistently-sized features.
#  -For a filter GFF, I'm using "-f 0.5" so stuff doesn't need to be fully
#   contained in the reference feature. This might be kind of wrong in some
#   cases like the exonic transcriptome. But it shouldn't be too bad.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <feature gff>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='filter_gff', help='Filter the TEs by overlap with genes in the given gff file [Default: %default]')
    parser.add_option('-r', dest='repeats_gff', default='%s/research/common/data/genomes/hg19/annotation/repeatmasker/hg19.fa.out.tp.gff' % os.environ['HOME'])
    parser.add_option('-s', dest='strand_split', default=False, action='store_true', help='Split statistics by strand [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide a gff file for the feature of interest.')
    else:
        feature_gff = args[0]

    ############################################
    # GFF filter
    ############################################
    # filter TEs and features by gff file
    if options.filter_gff:
        filter_merged_bed_fd, filter_merged_bed_file = tempfile.mkstemp()
        subprocess.call('sortBed -i %s | mergeBed -i - > %s' % (options.filter_gff, filter_merged_bed_file), shell=True)

        # filter TE GFF
        te_gff_fd, te_gff_file = tempfile.mkstemp()
        subprocess.call('intersectBed -a %s -b %s > %s' % (options.repeats_gff, filter_merged_bed_file, te_gff_file), shell=True)
        options.repeats_gff = te_gff_file

        # filter feature GFF
        feature_gff_gff_fd, feature_gff_gff_file = tempfile.mkstemp()
        subprocess.call('intersectBed -u -f 0.5 -a %s -b %s > %s' % (feature_gff, filter_merged_bed_file, feature_gff_gff_file), shell=True)
        feature_gff = feature_gff_gff_file

    ############################################
    # lengths
    ############################################
    # compute feature length
    feature_len, feature_num = feature_stats(feature_gff)    

    if feature_num == 0:
        print >> sys.stderr, 'Zero features'
        exit()

    # compute size of search space
    if options.filter_gff:
        genome_length = count_bed(filter_merged_bed_file, feature_len)
    else:
        genome_length = count_hg19()

    # hash counted repeat genomic bp
    te_lengths = te_target_size(options.repeats_gff, feature_len)

    ############################################
    # hash TE/feature overlaps
    ############################################
    # initialize
    te_features = {}
    for rep, fam in te_lengths:
        if options.strand_split:
            te_features[(rep+'+',fam)] = set()
            te_features[('*+',fam)] = set()
            te_features[('*+','*')] = set()
            te_features[(rep+'-',fam)] = set()
            te_features[('*-',fam)] = set()
            te_features[('*-','*')] = set()
        else:
            te_features[(rep,fam)] = set()
            te_features[('*',fam)] = set()
            te_features[('*','*')] = set()
        
    p = subprocess.Popen('intersectBed -wo -a %s -b %s' % (options.repeats_gff,feature_gff), shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        a = line.split('\t')
        
        kv = gff.gtf_kv(a[8])
        rep = kv['repeat']
        fam = kv['family']

        fchrom = a[9]
        fstart = int(a[12])
        fend = int(a[13])

        rep_star = '*'
        if options.strand_split:
            tstrand = a[6]
            fstrand = a[15]
            if tstrand == fstrand:
                rep += '+'
                rep_star += '+'
            else:
                rep += '-'
                rep_star += '-'

        te_features[(rep,fam)].add((fchrom,fstart,fend))
        te_features[(rep_star,fam)].add((fchrom,fstart,fend))
        te_features[(rep_star,'*')].add((fchrom,fstart,fend))

    p.communicate()

    ############################################SW
    # compute stats and print
    ############################################
    lines = []
    p_vals = []
    for te in te_features:
        rep, fam = te

        if options.strand_split:
            te_len = te_lengths[(rep[:-1],fam)]
            te_p = float(te_len) / (2*genome_length)
        else:
            te_len = te_lengths[(rep,fam)]
            te_p = float(te_len) / genome_length
        
        te_count = len(te_features.get(te,[]))
        exp_count = te_p * feature_num

        fold_change = te_count / exp_count

        if fold_change > 1:
            p_val = binom.sf(te_count-1, feature_num, te_p)
        else:
            p_val = binom.cdf(te_count, feature_num, te_p)
        
        p_vals.append(p_val)

        cols = (rep, fam, te_len, te_count, exp_count, fold_change, p_val)
        lines.append('%-18s %-18s %9d %8d %8.1f %8.2f %10.2e' % cols)

    # correct for multiple hypotheses correction
    q_vals = fdr.ben_hoch(p_vals)
    for i in range(len(lines)):
        qline = lines[i] + ' %10.2e' % q_vals[i]
        print qline

    ############################################
    # clean
    ############################################
    if options.filter_gff:
        os.close(filter_merged_bed_fd)
        os.remove(filter_merged_bed_file)
        os.close(te_gff_fd)
        os.remove(te_gff_file)
        os.close(feature_gff_gff_fd)
        os.remove(feature_gff_gff_file)


################################################################################
# feature_stats
#
# Compute mean feature length and count the number of features.
################################################################################
def feature_stats(feature_gff):
    feature_lengths = []
    for line in open(feature_gff):
        a = line.split('\t')
        feature_lengths.append(int(a[4])-int(a[3])+1)

    if len(feature_lengths) == 0:
        fmean = None
    else:
        fmean = int(0.5+stats.mean(feature_lengths))

    return fmean, len(feature_lengths)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
