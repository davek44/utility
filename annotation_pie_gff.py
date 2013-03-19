#!/usr/bin/env python
from optparse import OptionParser
import os, re, shutil, subprocess
import pysam
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.lib.ggplot2 as ggplot2

grdevices = importr('grDevices')

################################################################################
# annotation_pie_gff.py
#
# Count the # of GFF feature bp's overlapping various annotation classes and
# make pie charts.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <hg19|mm9> <gff>'
    parser = OptionParser(usage)
    parser.add_option('-t','--transcriptome', dest='transcriptome_only', default=False, action='store_true', help='Consider the transcriptome only, i.e. no intergenic class [Default: %default]')
    parser.add_option('-y', '--ylabel', dest='ylabel', default='Feature nt', help='Y-axis label [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) == 2:
        genome = args[0]
        gff_file = args[1]
    else:
        parser.error(usage)

    if genome == 'hg19':
        annotation_dir = '/n/rinn_data1/indexes/human/hg19/transcriptome/dk_pie'
        assembly_dir = '/n/home03/dkelley/research/common/data/genomes/hg19/assembly'
    elif genome == 'mm9':
        annotation_dir = '/n/rinn_data1/indexes/mouse/mm9/annotations/dk_pie'
        assembly_dir = '/n/home03/dkelley/research/common/data/genomes/mm9/assembly'
    else:
        parser.error('Genome must specify hg19 or mm9.')

    # count features
    gff_count = 0
    for line in open(gff_file):
        a = line.split('\t')
        gff_count += int(a[4]) - int(a[3]) + 1
    
    # CDS
    cds_length = annotation_length('%s/cds.bed' % annotation_dir, assembly_dir)
    cds_count = count_intersection(gff_file, '%s/cds.bed' % annotation_dir)
    cds_count_norm = cds_count / float(cds_length)

    # 3' UTR
    utr3_length = annotation_length('%s/utrs_3p.bed' % annotation_dir, assembly_dir)
    utr3_count = count_intersection(gff_file, '%s/utrs_3p.bed' % annotation_dir)
    utr3_count_norm = utr3_count / float(utr3_length)

    # 5' UTR
    utr5_length = annotation_length('%s/utrs_5p.bed' % annotation_dir, assembly_dir)
    utr5_count = count_intersection(gff_file, '%s/utrs_5p.bed' % annotation_dir)
    utr5_count_norm = utr5_count / float(utr5_length)

    # lncRNAs
    lnc_length = annotation_length('%s/lncrna.bed' % annotation_dir, assembly_dir)
    lnc_count = count_intersection(gff_file, '%s/lncrna.bed' % annotation_dir)
    lnc_count_norm = lnc_count / float(lnc_length)

    # intersect introns
    intron_length = annotation_length('%s/introns.bed' % annotation_dir, assembly_dir)
    intron_count = count_intersection(gff_file, '%s/introns.bed' % annotation_dir)
    intron_count_norm = intron_count / float(intron_length)

    # intergenic remains
    if not options.transcriptome_only:
        intergenic_length = genome_length(assembly_dir) - cds_length - utr3_length - utr5_length - lnc_length - intron_length
        intergenic_count = gff_count - cds_count - utr3_count - utr5_count - lnc_count - intron_count
        intergenic_count_norm = intergenic_count / float(intergenic_length)

    # print counts to file
    counts_out = open('counts.txt','w')
    print >> counts_out, 'CDS        %10d %8d %9.2e' % (cds_length, cds_count, cds_count_norm)
    print >> counts_out, "5'UTR      %10d %8d %9.2e" % (utr5_length, utr5_count, utr5_count_norm)
    print >> counts_out, "3'UTR      %10d %8d %9.2e" % (utr3_length, utr3_count, utr3_count_norm)
    print >> counts_out, 'lncRNA     %10d %8d %9.2e' % (lnc_length, lnc_count, lnc_count_norm)
    print >> counts_out, 'Intron     %10d %8d %9.2e' % (intron_length, intron_count, intron_count_norm)
    if not options.transcriptome_only:
        print >> counts_out, 'Intergenic %10d %8d %9.2e' % (intergenic_length, intergenic_count, intergenic_count_norm)
    counts_out.close()

    # construct data frame
    if options.transcriptome_only:
        dummy_r = ro.StrVector(['.']*5)
        counts_r = ro.IntVector([cds_count, utr3_count, utr5_count, lnc_count, intron_count])
        annotations_r = ro.StrVector(['CDS', '3\'UTR', '5\'UTR', 'lncRNA', 'Intron'])
    else:
        dummy_r = ro.StrVector(['.']*6)
        counts_r = ro.IntVector([cds_count, utr3_count, utr5_count, lnc_count, intron_count, intergenic_count])
        annotations_r = ro.StrVector(['CDS', '3\'UTR', '5\'UTR', 'lncRNA', 'Intron', 'Intergenic'])
    df = ro.DataFrame({'dummy':dummy_r, 'count':counts_r, 'annotation':annotations_r})

    # plot bar
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='dummy', y='count', fill='annotation') + \
        ggplot2.geom_bar(stat='identity', width=1) + \
        ggplot2.scale_x_discrete('') + \
        ggplot2.scale_y_continuous(options.ylabel) + \
        ggplot2.scale_fill_discrete('Annotation')

    # plot to file
    grdevices.pdf(file='annotation_gff_bar.pdf')
    gp.plot()
    grdevices.dev_off()

    # plot pie
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='dummy', y='count', fill='annotation') + \
        ggplot2.geom_bar(stat='identity', width=1) + \
        ggplot2.coord_polar(theta='y') + \
        ggplot2.scale_x_discrete('') + \
        ggplot2.scale_y_continuous(options.ylabel) + \
        ggplot2.scale_fill_discrete('Annotation')

    # plot to file
    grdevices.pdf(file='annotation_gff_pie.pdf')
    gp.plot()
    grdevices.dev_off()

    # normalize
    if options.transcriptome_only:
        counts_r = ro.FloatVector([cds_count_norm, utr3_count_norm, utr5_count_norm, lnc_count_norm, intron_count_norm])
    else:
        counts_r = ro.FloatVector([cds_count_norm, utr3_count_norm, utr5_count_norm, lnc_count_norm, intron_count_norm, intergenic_count_norm])
    df = ro.DataFrame({'dummy':dummy_r, 'count':counts_r, 'annotation':annotations_r})

    # plot norm bar
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='dummy', y='count', fill='annotation') + \
        ggplot2.geom_bar(stat='identity', width=1) + \
        ggplot2.scale_x_discrete('') + \
        ggplot2.scale_y_continuous('Length-normalized %s' % options.ylabel.lower()) + \
        ggplot2.scale_fill_discrete('Annotation')

    # plot to file
    grdevices.pdf(file='annotation_gff_norm_bar.pdf')
    gp.plot()
    grdevices.dev_off()

    # plot norm pie
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='dummy', y='count', fill='annotation') + \
        ggplot2.geom_bar(stat='identity', width=1) + \
        ggplot2.coord_polar(theta='y') + \
        ggplot2.scale_x_discrete('') + \
        ggplot2.scale_y_continuous('Length-normalized %s' % options.ylabel.lower()) + \
        ggplot2.scale_fill_discrete('Annotation')

    # plot to file
    grdevices.pdf(file='annotation_gff_norm_pie.pdf')
    gp.plot()
    grdevices.dev_off()


################################################################################
# annotation_length
#
# Input
#  bed_file: Annotation BED file
#
# Output
#  alength:  Summed lengths of the annotations.
################################################################################
def annotation_length(bed_file, assembly_dir):
    if assembly_dir.find('hg19') != -1:
        gaps_file = '%s/hg19_gaps.bed' % assembly_dir
    elif assembly_dir.find('mm9') != -1:
        gaps_file = '%s/mm9_gaps.bed' % assembly_dir
    else:
        print >> sys.stderr, 'Bad assembly directory'
        exit(1)

    alength = 0
    
    p = subprocess.Popen('subtractBed -a %s -b %s' % (bed_file,gaps_file), shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        a = line.split('\t')
        alength += int(a[2]) - int(a[1])
    p.communicate()
    
    return alength


################################################################################
# count_intersection
#
# Input
#  gff_file: Features GFF file
#  bed_file: Annotation BED file
#
# Output
#  count:    The number of nt overlapping the annotation in the BED file.
################################################################################
def count_intersection(gff_file, bed_file):
    # intersect
    p = subprocess.Popen('intersectBed -a %s -b %s' % (gff_file,bed_file), shell=True, stdout=subprocess.PIPE)

    # count
    count = 0
    for line in p.stdout:
        a = line.split('\t')
        count += int(a[4]) - int(a[3])

    return count

################################################################################
# genome_length
#
# Output
#  glength: Length of the human genome assembly minus gaps.
################################################################################
def genome_length(assembly_dir):
    if assembly_dir.find('hg19') != -1:
        chrom_file = '%s/human.hg19.genome' % assembly_dir
        gaps_file = '%s/hg19_gaps.bed' % assembly_dir
    elif assembly_dir.find('mm9') != -1:
        chrom_file = '%s/mouse.mm9.genome' % assembly_dir
        gaps_file = '%s/mm9_gaps.bed' % assembly_dir
    else:
        print >> sys.stderr, 'Bad assembly directory'
        exit(1)
        
    # count chromosome sizes
    glength = 0
    for line in open(chrom_file):
        a = line.split()
        chrom = a[0]
        if not chrom.startswith('chrUn') and chrom.find('random') == -1 and chrom.find('hap') == -1:
            glength += int(a[1])
            
    # subtract gaps
    for line in open(gaps_file):
        a = line.split()
        chrom = a[0]
        if not chrom.startswith('chrUn') and chrom.find('random') == -1 and chrom.find('hap') == -1:
            glength -= int(a[2]) - int(a[1])

    return glength


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
