#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import binom, norm
import gzip, os, subprocess, sys, tempfile
import fdr, gff, stats

################################################################################
# te_gff_ovl_enrich.py
#
# Compute statistics about the overlap of GFF/BED entries on TEs.
#
# Notes:
#  -My counting assumes that the GFF/BED features do not overlap.
################################################################################

home_dir = os.environ['HOME']

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <feature gff/bed>'
    parser = OptionParser(usage)
    parser.add_option('-r', dest='repeats_gff', default='%s/research/common/data/genomes/hg19/annotation/repeatmasker/hg19.fa.out.tp.gff' % home_dir)
    parser.add_option('-n', dest='null_iterations', type=int, default=50, help='Number of shuffles to perform to estimate null distribution [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide a gff file for the feature of interest.')
    else:
        feature_gff = args[0]        

    # count genomic bp
    genome_bp = count_hg19()

    # count feature bp
    if feature_gff[-3:] == 'bed':
        feature_bp = count_bed(feature_gff)
    else:
        feature_bp = count_gff(feature_gff)

    # hash counted repeat genomic bp
    te_in = open(options.repeats_gff)
    genome_te_bp = hash_te(te_in)
    te_in.close()

    # convert feature gff to bed
    if feature_gff[-3:] == 'gtf':
        feature_bed_fd, feature_bed_file = tempfile.mkstemp()
        subprocess.call('gtf2bed.py %s > %s' % (feature_gff,feature_bed_file), shell=True)
        

    elif feature_gff[-3:] == 'gff':
        feature_bed_fd, feature_bed_file = tempfile.mkstemp()
        subprocess.call('gff2bed.py %s > %s' % (feature_gff,feature_bed_file), shell=True)

    elif feature_gff[-3:] == 'bed':
        feature_bed_file = feature_gff
        
    else:
        parser.error('Cannot recognize gff format suffix')

    ############################################
    # null distribution
    ############################################
    shuffle_bed_fd, shuffle_bed_file = tempfile.mkstemp()

    te_null_bp = {}
    for ni in range(options.null_iterations):
        print >> sys.stderr, ni

        # shuffle feature bed
        p = subprocess.call('shuffleBed -i %s -g %s/research/common/data/genomes/hg19/assembly/human.hg19.genome -excl %s/research/common/data/genomes/hg19/assembly/hg19_gaps.bed > %s' % (feature_bed_file, home_dir, home_dir, shuffle_bed_file), shell=True)

        # intersect w/ TEs and hash overlaps
        te_tmp_bp = intersect_hash(options.repeats_gff, shuffle_bed_file)
        for te in genome_te_bp:
            te_null_bp.setdefault(te,[]).append(te_tmp_bp.get(te,0))

    ############################################
    # actual
    ############################################
    te_bp = intersect_hash(options.repeats_gff, feature_gff)

    ############################################
    # compute stats and print
    ############################################
    lines = []
    p_vals = []
    for te in genome_te_bp:
        feature_freq = float(te_bp.get(te,0))/feature_bp
        genome_freq = float(genome_te_bp[te])/genome_bp
        fold_change = feature_freq / genome_freq

        #print te, stats.mean(te_null_bp[te]), stats.sd(te_null_bp[te])

        null_u, null_sd = stats.mean_sd(te_null_bp[te])
        if null_sd == 0:
            null_sd = 1.0
            
        if fold_change > 1:
            p = norm.sf(te_bp[te]-1, loc=null_u, scale=null_sd)
        else:
            p = norm.cdf(te_bp.get(te,0), loc=null_u, scale=null_sd)

        p_vals.append(p)

        cols = (te[0], te[1], te_bp.get(te,0), feature_freq, genome_freq, fold_change, p)
        lines.append('%-18s %-18s %8d %11.2e %11.2e %9.2f %10.2e' % cols)

    # correct for multiple hypotheses correction
    q_vals = fdr.ben_hoch(p_vals)
    for i in range(len(lines)):
        qline = lines[i] + ' %10.2e' % q_vals[i]
        print qline

    ############################################
    # clean
    ############################################
    os.close(shuffle_bed_fd)
    os.remove(shuffle_bed_file)
    if feature_gff[-3:] != 'bed':
        os.close(feature_bed_fd)
        os.remove(feature_bed_file)


################################################################################
# count_hg19
#
# Count the number of bp in hg19 where TEs could be.
################################################################################
def count_hg19():
    chrom_sizes_file = '%s/research/common/data/genomes/hg19/assembly/human.hg19.genome' % home_dir
    gap_bed_file = '%s/research/common/data/genomes/hg19/assembly/hg19_gaps.bed' % home_dir
    valid_chrs = ['chr%d' % c for c in range(1,23)] + ['chrX','chrY']

    genome_bp = 0
    for line in open(chrom_sizes_file):        
        a = line.split()
        if len(a) > 0 and a[0] in valid_chrs:
            genome_bp += int(a[1])

    for line in open(gap_bed_file):
        a = line.split()
        if a[0] in valid_chrs:
            genome_bp -= int(a[2])-int(a[1])

    return genome_bp


################################################################################
# count_bed
#
# Count the number of covered nt in a bed file
################################################################################
def count_bed(bed_file):
    bp = 0
    for line in open(bed_file):
        a = line.split('\t')
        bp += int(a[2])-int(a[1])
    return bp


################################################################################
# count_gff
#
# Count the number of covered nt in a gff file
################################################################################
def count_gff(gff_file):
    bp = 0
    for line in open(gff_file):
        a = line.split('\t')
        bp += int(a[4])-int(a[3])+1
    return bp


################################################################################
# intersect_hash
#
# Intersect the TE GFF file with the feature GFF file and hash the overlap
# counts.
################################################################################
def intersect_hash(te_gff, feature_gff):
    te_bp = {}
    p = subprocess.Popen('intersectBed -a %s -b %s' % (te_gff, feature_gff), shell=True, stdout=subprocess.PIPE)
    te_bp = hash_te(p.stdout)
    p.communicate()
    return te_bp


################################################################################
# hash_te
#
# Hash the number of bp covered by various repeats in a RepeatMasker gff file or
# the intersection of a RepeatMasker gff file and the feature gtf file.
################################################################################
def hash_te(te_gff_in):
    te_bp = {}
    for line in te_gff_in:
        a = line.split('\t')

        kv = gff.gtf_kv(a[8])
        rep = kv['repeat']
        family = kv['family']

        length = int(a[4]) - int(a[3]) + 1

        te_bp[(rep,family)] = te_bp.get((rep,family),0) + length
        te_bp[('*',family)] = te_bp.get(('*',family),0) + length
        te_bp[('*','*')] = te_bp.get(('*','*'),0) + length
        if rep.startswith('LTR'):
            te_bp[('LTR*',family)] = te_bp.get(('LTR*',family),0) + length
        if rep.startswith('LTR12'):
            te_bp[('LTR12*',family)] = te_bp.get(('LTR12*',family),0) + length
        if rep.startswith('LTR7') and (len(rep) < 5 or rep[4].isalpha()):
            te_bp[('LTR7*',family)] = te_bp.get(('LTR7*',family),0) + length
        if rep.startswith('THE1') and len(rep) == 5:
            te_bp[('THE1*',family)] = te_bp.get(('THE1*',family),0) + length
        if rep.startswith('MER61') and len(rep) == 6:
            te_bp[('MER61*',family)] = te_bp.get(('MER61*',family),0) + length
        if rep.startswith('L1PA'):
            te_bp[('L1PA*',family)] = te_bp.get(('L1PA*',family),0) + length

    return te_bp


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
