#!/usr/bin/env python
from optparse import OptionParser
import os, re, shutil, subprocess
import pysam
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.lib.ggplot2 as ggplot2

grdevices = importr('grDevices')

################################################################################
# annotation_pie.py
#
# Count aligned reads to various annotation classes and make pie charts.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <hg19|mm9> <bam>'
    parser = OptionParser(usage)
    parser.add_option('-p', dest='paired', default=False, action='store_true', help='Paired end reads (only properly paired will be counted) [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) == 2:
        genome = args[0]
        bam_file = args[1]
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

    # count reads
    read_count = count_bam(bam_file, options.paired)
    
    # CDS
    cds_length = annotation_length('%s/cds.bed' % annotation_dir, assembly_dir)
    cds_reads = count_intersection(bam_file, '%s/cds.bed' % annotation_dir, options.paired)
    cds_reads_norm = cds_reads / float(cds_length)

    # 3' UTR
    utr3_length = annotation_length('%s/utrs_3p.bed' % annotation_dir, assembly_dir)
    utr3_reads = count_intersection(bam_file, '%s/utrs_3p.bed' % annotation_dir, options.paired)
    utr3_reads_norm = utr3_reads / float(utr3_length)

    # 5' UTR
    utr5_length = annotation_length('%s/utrs_5p.bed' % annotation_dir, assembly_dir)
    utr5_reads = count_intersection(bam_file, '%s/utrs_5p.bed' % annotation_dir, options.paired)
    utr5_reads_norm = utr5_reads / float(utr5_length)

    # lncRNAs
    lnc_length = annotation_length('%s/lncrna.bed' % annotation_dir, assembly_dir)
    lnc_reads = count_intersection(bam_file, '%s/lncrna.bed' % annotation_dir, options.paired)
    lnc_reads_norm = lnc_reads / float(lnc_length)

    # intersect introns
    intron_length = annotation_length('%s/introns.bed' % annotation_dir, assembly_dir)
    intron_reads = count_intersection(bam_file, '%s/introns.bed' % annotation_dir, options.paired)
    intron_reads_norm = intron_reads / float(intron_length)

    # intergenic remains
    intergenic_length = genome_length(assembly_dir) - cds_length - utr3_length - utr5_length - lnc_length - intron_length
    intergenic_reads = read_count - cds_reads - utr3_reads - utr5_reads - lnc_reads - intron_reads
    intergenic_reads_norm = intergenic_reads / float(intergenic_length)

    # print counts to file
    counts_out = open('counts.txt','w')
    print >> counts_out, 'CDS        %10d %8d %9.2e' % (cds_length, cds_reads, cds_reads_norm)
    print >> counts_out, "5'UTR      %10d %8d %9.2e" % (utr5_length, utr5_reads, utr5_reads_norm)
    print >> counts_out, "3'UTR      %10d %8d %9.2e" % (utr3_length, utr3_reads, utr3_reads_norm)
    print >> counts_out, 'lncRNA     %10d %8d %9.2e' % (lnc_length, lnc_reads, lnc_reads_norm)
    print >> counts_out, 'Intron     %10d %8d %9.2e' % (intron_length, intron_reads, intron_reads_norm)
    print >> counts_out, 'Intergenic %10d %8d %9.2e' % (intergenic_length, intergenic_reads, intergenic_reads_norm)
    counts_out.close()

    # construct data frame
    dummy_r = ro.StrVector(['.']*6)
    counts_r = ro.IntVector([cds_reads, utr3_reads, utr5_reads, lnc_reads, intron_reads, intergenic_reads])
    annotations_r = ro.StrVector(['CDS', '3\'UTR', '5\'UTR', 'lncRNA', 'Intron', 'Intergenic'])
    df = ro.DataFrame({'dummy':dummy_r, 'count':counts_r, 'annotation':annotations_r})

    # plot bar
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='dummy', y='count', fill='annotation') + \
        ggplot2.geom_bar(stat='identity', width=1) + \
        ggplot2.scale_x_discrete('') + \
        ggplot2.scale_y_continuous('Aligned reads') + \
        ggplot2.scale_fill_discrete('Annotation')

    # plot to file
    grdevices.pdf(file='annotation_reads_bar.pdf')
    gp.plot()
    grdevices.dev_off()

    # plot pie
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='dummy', y='count', fill='annotation') + \
        ggplot2.geom_bar(stat='identity', width=1) + \
        ggplot2.coord_polar(theta='y') + \
        ggplot2.scale_x_discrete('') + \
        ggplot2.scale_y_continuous('Aligned reads') + \
        ggplot2.scale_fill_discrete('Annotation')

    # plot to file
    grdevices.pdf(file='annotation_reads_pie.pdf')
    gp.plot()
    grdevices.dev_off()

    # normalize
    counts_r = ro.FloatVector([cds_reads_norm, utr3_reads_norm, utr5_reads_norm, lnc_reads_norm, intron_reads_norm, intergenic_reads_norm])
    df = ro.DataFrame({'dummy':dummy_r, 'count':counts_r, 'annotation':annotations_r})

    # plot norm bar
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='dummy', y='count', fill='annotation') + \
        ggplot2.geom_bar(stat='identity', width=1) + \
        ggplot2.scale_x_discrete('') + \
        ggplot2.scale_y_continuous('Length-normalized aligned reads') + \
        ggplot2.scale_fill_discrete('Annotation')

    # plot to file
    grdevices.pdf(file='annotation_reads_norm_bar.pdf')
    gp.plot()
    grdevices.dev_off()

    # plot norm pie
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='dummy', y='count', fill='annotation') + \
        ggplot2.geom_bar(stat='identity', width=1) + \
        ggplot2.coord_polar(theta='y') + \
        ggplot2.scale_x_discrete('') + \
        ggplot2.scale_y_continuous('Length-normalized aligned reads') + \
        ggplot2.scale_fill_discrete('Annotation')

    # plot to file
    grdevices.pdf(file='annotation_reads_norm_pie.pdf')
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
# count_bam
#
# Input
#  bam_file:
#  paired:     Reads are paired end, in which case count only properly paired.
#
# Output
#  read_count: Number of aligned fragments.
################################################################################
def count_bam(bam_file, paired):
    read_count = 0
    bam_in = pysam.Samfile(bam_file, 'rb')
    for aligned_read in bam_in:
        chrom = bam_in.getrname(aligned_read.tid)
        if not chrom.startswith('chrUn') and chrom.find('random') == -1 and chrom.find('hap') == -1:
            if not paired:
                read_count += 1.0/aligned_read.opt('NH')
            else:
                if aligned_read.is_proper_pair:
                    read_count += 0.5/aligned_read.opt('NH')
    return read_count


################################################################################
# count_intersection
#
# Input
#  bam_file: Read alignment BAM file
#  bed_file: Annotation BED file
#  paired:   Reads are paired end, in which case count only properly paired.
#
# Output
#  reads:    The number of reads (corrected for multi-mappers) overlapping the
#             annotation in the BED file.
################################################################################
def count_intersection(bam_file, bed_file, paired):
    # intersect
    subprocess.call('intersectBed -abam %s -b %s > tmp.bam' % (bam_file,bed_file), shell=True)

    # count
    reads = count_bam('tmp.bam', paired)

    # clean
    os.remove('tmp.bam')
    
    return reads
    

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
