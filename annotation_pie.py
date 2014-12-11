#!/usr/bin/env python
from optparse import OptionParser
import math, os, re, sys, subprocess, tempfile
import pysam
import ggplot

################################################################################
# annotation_bars.py
#
# Count aligned reads to various annotation classes and make pie charts.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <hg19|mm9> <bam1,bam2,...>'
    parser = OptionParser(usage)
    parser.add_option('-a', dest='annotations', default='rrna,smallrna,cds,utrs_3p,utrs_5p,pseudogene,lncrna,introns,intergenic', help='Comma-separated list of annotation classes to include [Default: %default]')
    parser.add_option('-o', dest='output_prefix', default='annotation', help='Output file prefix [Default: %default]')
    parser.add_option('-p', dest='paired_stranded', action='store_true', default=False, help='Paired end stranded reads, so split intersects by XS tag and strand [Default: %default]')
    parser.add_option('-t', dest='title', default='', help='Plot title [Default: %default]')
    parser.add_option('-u', dest='unstranded', action='store_true', default=False, help='Unstranded reads, so count intergenic and renormalize to lessen the impact of double counting [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) == 2:
        genome = args[0]
        bam_files = args[1].split(',')
    else:
        parser.error(usage)

    if genome == 'hg19':
        assembly_dir = '%s/research/common/data/genomes/hg19/assembly' % os.environ['HOME']
        if options.paired_stranded:
            annotation_dir = '%s/pie_stranded' % os.environ['GENCODE']
        else:
            annotation_dir = '%s/pie_unstranded' % os.environ['GENCODE']

    elif genome == 'mm9':
        assembly_dir = '%s/research/common/data/genomes/mm9/assembly' % os.environ['HOME']
        if options.paired_stranded:
            print >> sys.stderr, 'Stranded annotation BEDs dont exist for mm9'
            exit(1)
        else:
            annotation_dir = '/n/rinn_data1/indexes/mouse/mm9/annotations/dk_pie'

    else:
        parser.error('Genome must specify hg19 or mm9.')

    if options.paired_stranded:
        # split bam files by strand
        for bam_file in bam_files:
            split_bam_xs(bam_file)

    annotation_classes = set(options.annotations.split(','))

    ############################################
    # annotation lengths
    ############################################
    genome_length = count_genome(assembly_dir)

    annotation_lengths = {}
    for ann in annotation_classes:
        if ann != 'intergenic':
            annotation_bed = '%s/%s.bed' % (annotation_dir,ann)
            if os.path.isfile(annotation_bed):
                annotation_lengths[ann] = annotation_length(annotation_bed, assembly_dir)
            else:
                parser.error('Cannot find annotation BED %s' % annotation_bed)
                
    if 'intergenic' in annotation_classes:
        other_annotations_summed = sum(annotation_lengths.values())
        annotation_lengths['intergenic'] = genome_length - other_annotations_summed

    ############################################
    # annotation read counts
    ############################################
    genome_reads = 0
    for bam_file in bam_files:
        genome_reads += count_bam(bam_file)

    annotation_reads = {}
    for ann in annotation_classes:
        if ann != 'intergenic':
            annotation_bed = '%s/%s.bed' % (annotation_dir,ann)
            annotation_reads[ann] = 0
            for bam_file in bam_files:
                annotation_reads[ann] += count_intersection(bam_file, annotation_bed, options.unstranded, options.paired_stranded)

    if 'intergenic' in annotation_classes:
        other_annotations_summed = sum(annotation_reads.values())
        annotation_reads['intergenic'] = genome_reads - other_annotations_summed
    
        if options.unstranded:
            intergenic_reads_sub = annotation_reads['intergenic']
            intergenic_reads = 0
            for bam_file in bam_files:
                intergenic_reads += count_sans_intersection(bam_file, '%s/../gencode.v18.annotation.prerna.gtf' % annotation_dir)

    if options.paired_stranded:
        for bam_file in bam_files:
            os.remove(bam_file[:-4] + '_p.bam')
            os.remove(bam_file[:-4] + '_m.bam')

    ############################################
    # table
    ############################################
    annotation_labels = {'rrna':'rRNA', 'smallrna':'smallRNA', 'cds':'CDS', 'utrs_3p':'3\'UTR', 'utrs_5p':'5\'UTR', 'pseudogene':'Pseudogene', 'lncrna':'lncRNA', 'exons':'Exons', 'introns':'Introns', 'intergenic':'Intergenic', 'mrna':'mRNA'}

    reads_sum = float(sum(annotation_reads.values()))
    lengths_sum = float(sum(annotation_lengths.values()))

    annotation_ratio = {}

    counts_out = open('%s_counts.txt' % options.output_prefix, 'w')
    for ann in annotation_classes:
        read_pct = annotation_reads[ann]/reads_sum
        length_pct = annotation_lengths[ann]/lengths_sum

        if read_pct > 0:
            annotation_ratio[ann] = math.log(read_pct/length_pct,2)
        else:
            annotation_ratio[ann] = math.log((1+annotation_reads[ann])/(1+reads_sum),2)

        cols = (annotation_labels[ann], annotation_reads[ann], read_pct, length_pct, annotation_ratio[ann])
        print >> counts_out, '%10s  %8d  %6.4f  %6.4f  %5.2f' % cols
    counts_out.close()

    ############################################
    # pie chart
    ############################################
    pie_df = {'dummy':[], 'annotation':[], 'count':[]}
    for ann in annotation_classes:
        pie_df['dummy'].append('.')
        pie_df['annotation'].append(annotation_labels[ann])
        pie_df['count'].append(annotation_reads[ann])

    out_pdf = '%s_pie.pdf' % options.output_prefix
    ggplot.plot('%s/annotation_pie_pie.r'%os.environ['RDIR'], pie_df, [options.title, out_pdf], df_file=out_pdf[:-1])

    ############################################
    # read:length ratio
    ############################################
    ratio_df = {'annotation':[], 'ratio':[]}
    for ann in annotation_classes:
        ratio_df['annotation'].append(annotation_labels[ann])
        ratio_df['ratio'].append(annotation_ratio[ann])

    ggplot.plot('%s/annotation_pie_ratios.r'%os.environ['RDIR'], ratio_df, [options.title, '%s_ratios.pdf'%options.output_prefix])


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
#  bam_file:        Read alignment BAM file
#  bed_file:        Annotation BED file
#  paired_stranded: Reads are paired and stranded
#  introns:         Counting introns, so ignore spliced reads
#
# Output
#  reads:           The number of reads (corrected for multi-mappers)
#                    overlapping the annotation in the BED file.
################################################################################
def count_intersection(bam_file, bed_file, unstranded, paired_stranded, introns=False):
    intersect_bam_fd, intersect_bam_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])

    if paired_stranded:
        # make split bed temp files
        bedp_fd, bedp_file = tempfile.mkstemp()
        bedm_fd, bedm_file = tempfile.mkstemp()

        # split bed file
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
# count_genome
#
# Output
#  glength: Length of the human genome assembly minus gaps.
################################################################################
def count_genome(assembly_dir):
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
