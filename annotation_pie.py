#!/usr/bin/env python
from optparse import OptionParser
import os, re, shutil, subprocess, tempfile
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
    parser.add_option('-p', dest='paired', action='store_true', default=False, help='Paired end reads, so split intersects by XS tag and strand [Default: %default]')
    parser.add_option('-u', dest='unstranded', action='store_true', default=False, help='Unstranded reads, so count intergenic and renormalize to lessen the impact of double counting [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) == 2:
        genome = args[0]
        bam_file = args[1]
    else:
        parser.error(usage)

    if genome == 'hg19':
        annotation_dir = '%s/research/common/data/genomes/hg19/annotation/gencode_v15/pie' % os.environ['HOME']
        assembly_dir = '%s/research/common/data/genomes/hg19/assembly' % os.environ['HOME']
    elif genome == 'mm9':
        annotation_dir = '/n/rinn_data1/indexes/mouse/mm9/annotations/dk_pie'
        assembly_dir = '%s/research/common/data/genomes/mm9/assembly' % os.environ['HOME']
    else:
        parser.error('Genome must specify hg19 or mm9.')

    # count reads
    read_count = count_bam(bam_file)

    if options.paired:
        # split bam file by strand
        split_bam_xs(bam_file)
    
    # rRNA
    rrna_length = annotation_length('%s/rrna.bed' % annotation_dir, assembly_dir)
    rrna_reads = count_intersection(bam_file, '%s/rrna.bed' % annotation_dir, options.unstranded, options.paired)
    rrna_reads_norm = rrna_reads / float(rrna_length)

    # small RNA
    smallrna_length = annotation_length('%s/smallrna.bed' % annotation_dir, assembly_dir)
    smallrna_reads = count_intersection(bam_file, '%s/smallrna.bed' % annotation_dir, options.unstranded, options.paired)
    smallrna_reads_norm = smallrna_reads / float(smallrna_length)

    # CDS
    cds_length = annotation_length('%s/cds.bed' % annotation_dir, assembly_dir)
    cds_reads = count_intersection(bam_file, '%s/cds.bed' % annotation_dir, options.unstranded, options.paired)
    cds_reads_norm = cds_reads / float(cds_length)

    # 3' UTR
    utr3_length = annotation_length('%s/utrs_3p.bed' % annotation_dir, assembly_dir)
    utr3_reads = count_intersection(bam_file, '%s/utrs_3p.bed' % annotation_dir, options.unstranded, options.paired)
    utr3_reads_norm = utr3_reads / float(utr3_length)

    # 5' UTR
    utr5_length = annotation_length('%s/utrs_5p.bed' % annotation_dir, assembly_dir)
    utr5_reads = count_intersection(bam_file, '%s/utrs_5p.bed' % annotation_dir, options.unstranded, options.paired)
    utr5_reads_norm = utr5_reads / float(utr5_length)

    # pseudogenes
    pseudo_length = annotation_length('%s/pseudogene.bed' % annotation_dir, assembly_dir)
    pseudo_reads = count_intersection(bam_file, '%s/pseudogene.bed' % annotation_dir, options.unstranded, options.paired)
    pseudo_reads_norm = pseudo_reads / float(pseudo_length)

    # lncRNAs
    lnc_length = annotation_length('%s/lncrna.bed' % annotation_dir, assembly_dir)
    lnc_reads = count_intersection(bam_file, '%s/lncrna.bed' % annotation_dir, options.unstranded, options.paired)
    lnc_reads_norm = lnc_reads / float(lnc_length)

    # intersect introns
    intron_length = annotation_length('%s/introns.bed' % annotation_dir, assembly_dir)
    intron_reads = count_intersection(bam_file, '%s/introns.bed' % annotation_dir, options.unstranded, options.paired, introns=True)
    intron_reads_norm = intron_reads / float(intron_length)

    # intergenic
    intergenic_length = genome_length(assembly_dir) - rrna_length - smallrna_length - cds_length - utr3_length - utr5_length - pseudo_length - lnc_length - intron_length
    intergenic_reads = read_count - rrna_reads - smallrna_reads - cds_reads - utr3_reads - utr5_reads - pseudo_reads - lnc_reads - intron_reads
    if options.unstranded:
        intergenic_reads_sub = intergenic_reads
        intergenic_reads = count_sans_intersection(bam_file, '%s/../gencode.v15.annotation.gtf' % annotation_dir)

        read_count = rrna_reads + smallrna_reads + cds_reads + utr3_reads + utr5_reads + pseudo_reads + lnc_reads + intron_reads + intergenic_reads
    intergenic_reads_norm = intergenic_reads / float(intergenic_length)

    if options.paired:
        os.remove(bam_file[:-4] + '_p.bam')
        os.remove(bam_file[:-4] + '_m.bam')

    # compute pie values to print
    rrna_pie = float(rrna_reads)/read_count
    smallrna_pie = float(smallrna_reads)/read_count
    cds_pie = float(cds_reads)/read_count
    utr5_pie = float(utr5_reads)/read_count
    utr3_pie = float(utr3_reads)/read_count
    pseudo_pie = float(pseudo_reads)/read_count
    lnc_pie = float(lnc_reads)/read_count
    intron_pie = float(intron_reads)/read_count
    intergenic_pie = float(intergenic_reads)/read_count

    norm_sum = rrna_reads_norm + smallrna_reads_norm + cds_reads_norm + utr5_reads_norm + utr3_reads_norm + pseudo_reads_norm + lnc_reads_norm + intron_reads_norm + intergenic_reads_norm
    rrna_norm_pie = float(rrna_reads_norm)/norm_sum
    smallrna_norm_pie = float(smallrna_reads_norm)/norm_sum
    cds_norm_pie = float(cds_reads_norm)/norm_sum
    utr5_norm_pie = float(utr5_reads_norm)/norm_sum
    utr3_norm_pie = float(utr3_reads_norm)/norm_sum
    pseudo_norm_pie = float(pseudo_reads_norm)/norm_sum
    lnc_norm_pie = float(lnc_reads_norm)/norm_sum
    intron_norm_pie = float(intron_reads_norm)/norm_sum
    intergenic_norm_pie = float(intergenic_reads_norm)/norm_sum

    # print counts to file
    counts_out = open('annotation_counts.txt','w')
    print >> counts_out, 'rRNA        %10d %8d %5.3f %9.2e %5.3f' % (rrna_length, rrna_reads, rrna_pie, rrna_reads_norm, rrna_norm_pie)
    print >> counts_out, 'smallRNA    %10d %8d %5.3f %9.2e %5.3f' % (smallrna_length, smallrna_reads, smallrna_pie, smallrna_reads_norm, smallrna_norm_pie)
    print >> counts_out, 'CDS         %10d %8d %5.3f %9.2e %5.3f' % (cds_length, cds_reads, cds_pie, cds_reads_norm, cds_norm_pie)
    print >> counts_out, "5'UTR       %10d %8d %5.3f %9.2e %5.3f" % (utr5_length, utr5_reads, utr5_pie, utr5_reads_norm, utr5_norm_pie)
    print >> counts_out, "3'UTR       %10d %8d %5.3f %9.2e %5.3f" % (utr3_length, utr3_reads, utr3_pie, utr3_reads_norm, utr3_norm_pie)
    print >> counts_out, 'Pseudogene  %10d %8d %5.3f %9.2e %5.3f' % (pseudo_length, pseudo_reads, pseudo_pie, pseudo_reads_norm, pseudo_norm_pie)
    print >> counts_out, 'lncRNA      %10d %8d %5.3f %9.2e %5.3f' % (lnc_length, lnc_reads, lnc_pie, lnc_reads_norm, lnc_norm_pie)
    print >> counts_out, 'Intron      %10d %8d %5.3f %9.2e %5.3f' % (intron_length, intron_reads, intron_pie, intron_reads_norm, intron_norm_pie)
    print >> counts_out, 'Intergenic  %10d %8d %5.3f %9.2e %5.3f' % (intergenic_length, intergenic_reads, intergenic_pie, intergenic_reads_norm, intergenic_norm_pie)
    if options.unstranded:
        print >> counts_out, 'IntergenSub %10d %8d' % (intergenic_length, intergenic_reads_sub)
    counts_out.close()

    # construct data frame
    dummy_r = ro.StrVector(['.']*9)
    counts_r = ro.IntVector([rrna_reads, smallrna_reads, cds_reads, utr3_reads, utr5_reads, pseudo_reads, lnc_reads, intron_reads, intergenic_reads])
    annotations_r = ro.StrVector(['rRNA', 'smallRNA', 'CDS', '3\'UTR', '5\'UTR', 'lncRNA', 'Pseudogene', 'Intron', 'Intergenic'])
    df = ro.DataFrame({'dummy':dummy_r, 'count':counts_r, 'annotation':annotations_r})

    # plot bar
    '''
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
    '''

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
    counts_r = ro.FloatVector([rrna_reads_norm, smallrna_reads_norm, cds_reads_norm, utr3_reads_norm, utr5_reads_norm, pseudo_reads_norm, lnc_reads_norm, intron_reads_norm, intergenic_reads_norm])
    df = ro.DataFrame({'dummy':dummy_r, 'count':counts_r, 'annotation':annotations_r})

    # plot norm bar
    '''
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
    '''

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

    # normalize, sans small RNA
    norm_sum = cds_reads_norm + utr5_reads_norm + utr3_reads_norm + pseudo_reads_norm + lnc_reads_norm + intron_reads_norm + intergenic_reads_norm
    cds_norm_pie = float(cds_reads_norm)/norm_sum
    utr5_norm_pie = float(utr5_reads_norm)/norm_sum
    utr3_norm_pie = float(utr3_reads_norm)/norm_sum
    pseudo_norm_pie = float(pseudo_reads_norm)/norm_sum
    lnc_norm_pie = float(lnc_reads_norm)/norm_sum
    intron_norm_pie = float(intron_reads_norm)/norm_sum
    intergenic_norm_pie = float(intergenic_reads_norm)/norm_sum

    counts_r = ro.FloatVector([cds_reads_norm, utr3_reads_norm, utr5_reads_norm, pseudo_reads_norm, lnc_reads_norm, intron_reads_norm, intergenic_reads_norm])
    annotations_r = ro.StrVector(['CDS', '3\'UTR', '5\'UTR', 'lncRNA', 'Pseudogene', 'Intron', 'Intergenic'])
    df = ro.DataFrame({'dummy':dummy_r, 'count':counts_r, 'annotation':annotations_r})

    # plot norm bar
    '''
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='dummy', y='count', fill='annotation') + \
        ggplot2.geom_bar(stat='identity', width=1) + \
        ggplot2.scale_x_discrete('') + \
        ggplot2.scale_y_continuous('Length-normalized aligned reads') + \
        ggplot2.scale_fill_discrete('Annotation')

    # plot to file
    grdevices.pdf(file='annotation_reads_normr_bar.pdf')
    gp.plot()
    grdevices.dev_off()
    '''

    # plot norm pie
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='dummy', y='count', fill='annotation') + \
        ggplot2.geom_bar(stat='identity', width=1) + \
        ggplot2.coord_polar(theta='y') + \
        ggplot2.scale_x_discrete('') + \
        ggplot2.scale_y_continuous('Length-normalized aligned reads') + \
        ggplot2.scale_fill_discrete('Annotation')

    # plot to file
    grdevices.pdf(file='annotation_reads_normr_pie.pdf')
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
#
# Output
#  read_count: Number of aligned fragments.
################################################################################
def count_bam(bam_file, skip_spliced=False):
    read_count = 0

    bam_in = pysam.Samfile(bam_file, 'rb')
    for aligned_read in bam_in:
        # high quality
        if aligned_read.mapq > 0:
            # we're not skipping spliced or it's not spliced
            if not skip_spliced or not spliced(aligned_read):
                # not a dumb chromosome
                chrom = bam_in.getrname(aligned_read.tid)
                if not chrom.startswith('chrUn') and chrom.find('random') == -1 and chrom.find('hap') == -1:
                    if aligned_read.is_paired:
                        read_count += 0.5/aligned_read.opt('NH')
                    else:
                        read_count += 1.0/aligned_read.opt('NH')
    bam_in.close()

    return read_count


################################################################################
# count_intersection
#
# Input
#  bam_file: Read alignment BAM file
#  bed_file: Annotation BED file
#  paired:   Reads are paired
#  introns:  Counting introns, so ignore spliced reads
#
# Output
#  reads:    The number of reads (corrected for multi-mappers) overlapping the
#             annotation in the BED file.
################################################################################
def count_intersection(bam_file, bed_file, unstranded, paired, introns=False):
    intersect_bam_fd, intersect_bam_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])

    if paired:
        # make split bed temp files
        bedp_fd, bedp_file = tempfile.mkstemp()
        bedm_fd, bedm_file = tempfile.mkstemp()

        # splite bed file
        split_bed(bed_file, bedp_file, bedm_file)

        # get split BAM file names
        bamp_file = bam_file[:-4] + '_p.bam'
        bamm_file = bam_file[:-4] + '_m.bam'        

        # count +
        subprocess.call('intersectBed -f 0.5 -abam %s -b %s > %s' % (bamp_file,bedp_file,intersect_bam_file), shell=True)
        readsp = count_bam(intersect_bam_file, skip_spliced=introns)

        # count -
        subprocess.call('intersectBed -f 0.5 -abam %s -b %s > %s' % (bamm_file,bedm_file,intersect_bam_file), shell=True)
        readsm = count_bam(intersect_bam_file, skip_spliced=introns)

        # sum + and -
        reads = readsp + readsm

        # clean temp
        os.close(bedp_fd)
        os.remove(bedp_file)
        os.close(bedm_fd)
        os.remove(bedm_file)
    else:
        strand_str = '-s'
        if unstranded:
            strand_str = ''

        # intersect
        subprocess.call('intersectBed -f 0.5 %s -abam %s -b %s > %s' % (strand_str,bam_file,bed_file,intersect_bam_file), shell=True)

        # count
        reads = count_bam(intersect_bam_file, skip_spliced=introns)

    # clean
    os.close(intersect_bam_fd)
    os.remove(intersect_bam_file)
    
    return reads


################################################################################
# count_sans_intersection
#
# Input
#  bam_file: Read alignment BAM file
#  bed_file: Annotation BED file
#
# Output
#  reads:    The number of reads (corrected for multi-mappers) not overlapping
#            the annotation in the BED file.
################################################################################
def count_sans_intersection(bam_file, bed_file):
    intersect_bam_fd, intersect_bam_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])

    # intersect
    subprocess.call('intersectBed -s -v -f 0.5 -abam %s -b %s > %s' % (bam_file,bed_file,intersect_bam_file), shell=True)

    # count
    reads = count_bam(intersect_bam_file)

    # clean
    os.close(intersect_bam_fd)
    os.remove(intersect_bam_file)
    
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
# spliced
#
# Return true if the read is spliced.
################################################################################
def spliced(aligned_read):
    spliced = False
    for code,size in aligned_read.cigar:
        if code == 3:
            spliced = True
    return spliced


################################################################################
# split_bam_xs
#
# Split the alignments in the given BAM file into plus and minus strand.
################################################################################
def split_bam_xs(bam_file):
    bamp_file = bam_file[:-4] + '_p.bam'
    bamm_file = bam_file[:-4] + '_m.bam'

    bam_in = pysam.Samfile(bam_file, 'rb')
    bamp_out = pysam.Samfile(bamp_file, 'wb', template=bam_in)
    bamm_out = pysam.Samfile(bamm_file, 'wb', template=bam_in)
    
    for aligned_read in bam_in:
        if aligned_read.opt('XS') == '+':
            bamp_out.write(aligned_read)
        else:
            bamm_out.write(aligned_read)

    bam_in.close()
    bamp_out.close()
    bamm_out.close()


################################################################################
# split_bed
################################################################################
def split_bed(bed_file, bedp_file, bedm_file):
    bedp_out = open(bedp_file, 'w')
    bedm_out = open(bedm_file, 'w')

    for line in open(bed_file):
        a = line.split('\t')
        if a[5] == '+':
            print >> bedp_out, line,
        else:
            print >> bedm_out, line,

    bedp_out.close()
    bedm_out.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
