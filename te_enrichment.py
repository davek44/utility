#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import binom, norm
import gzip, os, subprocess, sys
import fdr, gff, stats

################################################################################
# te_enrichment.py
#
# Compute a p-value for the enrichment of each TE family in a given GFF/BED file.
################################################################################

null_iterations = 100
home_dir = os.environ['HOME']

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <feature gff/bed>'
    parser = OptionParser(usage)
    parser.add_option('-f', dest='feature', help='Use only rows from the feature gff where the 3rd column matches feature')
    parser.add_option('-r', dest='repeats_gff', default='%s/research/common/data/genomes/hg19/annotation/repeatmasker/hg19.fa.out.tp.gff' % home_dir)
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide a gff file for the feature of interest.')
    else:
        feature_gff = args[0]        

    # filter for exons
    filt_feature_gff = feature_gff
    if options.feature:
        filt_feature_gff += '.tmp'
        subprocess.call("awk '$3 == \"%s\"' %s > %s" % (options.feature,feature_gff,filt_feature_gff), shell=True)

    # count genomic bp
    genome_bp = count_hg19()

    # count feature bp
    if filt_feature_gff[-3:] == 'bed':
        feature_bp = count_bed(filt_feature_gff)
    else:
        feature_bp = count_gff(filt_feature_gff)

    # hash counted repeat genomic bp
    genome_te_bp = hash_te(options.repeats_gff)

    # convert feature gff to bed
    feature_pre = os.path.splitext(feature_gff)[0]
    if feature_gff[-3:] == 'gtf':
        if options.feature == 'CDS':
            subprocess.call('gtf2bed.py -c %s > %s.bed' % (filt_feature_gff,feature_pre), shell=True)
        else:
            subprocess.call('gtf2bed.py %s > %s.bed' % (filt_feature_gff,feature_pre), shell=True)

    elif feature_gff[-3:] == 'gff':
        subprocess.call('gff2bed.py %s > %s.bed' % (filt_feature_gff,feature_pre), shell=True)
        
    elif feature_gff[-3:] != 'bed':
        parser.error('Cannot recognize gff format suffix')

    te_null_bp = {}
    for i in range(null_iterations):
        print >> sys.stderr, i

        # shuffle feature bed
        p = subprocess.call('shuffleBed -i %s.bed -g %s/research/common/data/genomes/hg19/assembly/human.hg19.genome -excl %s/research/common/data/genomes/hg19/assembly/hg19_gaps.bed > shuffle.bed' % (feature_pre, home_dir, home_dir), shell=True)

        # intersect w/ TEs
        subprocess.call('intersectBed -a %s -b shuffle.bed > te_shuffle.gff' % options.repeats_gff, shell=True)

        # hash by TE
        te_tmp_bp = hash_te('te_shuffle.gff')
        for te in genome_te_bp:
            te_null_bp.setdefault(te,[]).append(te_tmp_bp.get(te,0))

    # intersect real w/ TEs
    subprocess.call('intersectBed -a %s -b %s > te_feature.gff' % (options.repeats_gff, filt_feature_gff), shell=True)

    # hash by TE
    te_bp = hash_te('te_feature.gff')

    if not os.path.isdir('null_out'):
        os.mkdir('null_out')

    # print statistics
    lines = []
    p_vals = []
    for te in genome_te_bp:
        feature_freq = float(te_bp.get(te,0))/feature_bp
        genome_freq = float(genome_te_bp[te])/genome_bp
        fold_change = feature_freq / genome_freq

        null_out = open('null_out/%s_%s.txt' % (te[0].replace('/','_'), te[1].replace('/','_')), 'w')
        print >> null_out, '\n'.join([str(nb) for nb in te_null_bp[te]])
        null_out.close()

        #print te, stats.mean(te_null_bp[te]), stats.sd(te_null_bp[te])

        if fold_change > 1:
            p = norm.sf(te_bp[te]-1, loc=stats.mean(te_null_bp[te]), scale=stats.sd(te_null_bp[te]))
            #p = len([null_bp for null_bp in te_null_bp[te] if te_bp[te] <= null_bp])/float(null_iterations)
        else:
            p = norm.cdf(te_bp.get(te,0), loc=stats.mean(te_null_bp[te]), scale=stats.sd(te_null_bp[te]))
            #p = len([null_bp for null_bp in te_null_bp[te] if te_bp.get(te,0) >= null_bp])/float(null_iterations)

        if stats.mean(te_null_bp[te]) > 0:
            p_vals.append(p)

            cols = (te[0], te[1], te_bp.get(te,0), feature_freq, genome_freq, fold_change, p)
            lines.append('%-18s %-18s %8d %11.2e %11.2e %9.2f %10.2e' % cols)

    # correct for multiple hypotheses correction
    q_vals = fdr.ben_hoch(p_vals)
    for i in range(len(lines)):
        qline = lines[i] + ' %10.2e' % q_vals[i]
        print qline


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
# hash_te
#
# Hash the number of bp covered by various repeats in a RepeatMasker gff file or
# the intersection of a RepeatMasker gff file and the feature gtf file.
################################################################################
def hash_te(te_gff_file):
    te_bp = {}
    for line in open(te_gff_file):
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
