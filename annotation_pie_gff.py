#!/usr/bin/env python
from optparse import OptionParser
import math, os, re, shutil, subprocess
import pysam
import annotation_pie, ggplot

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
    parser.add_option('-a', dest='annotations', default='cds,utrs_3p,utrs_5p,lncrna,introns', help='Comma-separated list of annotation classes to include [Default: %default]')
    parser.add_option('-o', dest='output_prefix', default='annotation', help='Output file prefix [Default: %default]')
    parser.add_option('-t', dest='title', default='Title', help='Plot title [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) == 2:
        genome = args[0]
        gff_file = args[1]
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

    annotation_classes = set(options.annotations.split(','))

    ############################################
    # annotation lengths
    ############################################
    genome_length = annotation_pie.count_genome(assembly_dir)

    annotation_lengths = {}
    for ann in annotation_classes:
        if ann != 'intergenic':
            annotation_bed = '%s/%s.bed' % (annotation_dir,ann)
            if os.path.isfile(annotation_bed):
                annotation_lengths[ann] = annotation_pie.annotation_length(annotation_bed, assembly_dir)
            else:
                parser.error('Cannot find annotation BED %s' % annotation_bed)
                
    if 'intergenic' in annotation_classes:
        other_annotations_summed = sum(annotation_lengths.values())
        annotation_lengths['intergenic'] = genome_length - other_annotations_summed

    ############################################
    # annotation feature counts
    ############################################
    genome_features = int(subprocess.check_output('wc -l %s' % gff_file, shell=True).split()[0])

    annotation_features = {}
    for ann in annotation_classes:
        if ann != 'intergenic':
            annotation_bed = '%s/%s.bed' % (annotation_dir,ann)
            annotation_features[ann] = count_intersection(gff_file, annotation_bed)

    if 'intergenic' in annotation_classes:
        other_annotations_summed = sum(annotation_features.values())
        annotation_features['intergenic'] = genome_reads - other_annotations_summed        
    
    ############################################
    # table
    ###########################################
    annotation_labels = {'rrna':'rRNA', 'smallrna':'smallRNA', 'cds':'CDS', 'utrs_3p':'3\'UTR', 'utrs_5p':'5\'UTR', 'pseudogene':'Pseudogene', 'lncrna':'lncRNA', 'introns':'Introns', 'intergenic':'Intergenic'}

    features_sum = float(sum(annotation_features.values()))
    lengths_sum = float(sum(annotation_lengths.values()))

    annotation_ratio = {}

    counts_out = open('%s_counts.txt' % options.output_prefix, 'w')
    for ann in annotation_classes:
        feature_pct = annotation_features[ann]/features_sum
        length_pct = annotation_lengths[ann]/lengths_sum

        if feature_pct > 0:
            annotation_ratio[ann] = math.log(feature_pct/length_pct,2)
        else:
            annotation_ratio[ann] = math.log((1+annotation_features[ann])/(1+features_sum),2)

        cols = (annotation_labels[ann], annotation_features[ann], feature_pct, length_pct, annotation_ratio[ann])
        print >> counts_out, '%10s  %8d  %6.4f  %6.4f  %5.2f' % cols
    counts_out.close()

    ############################################
    # pie chart
    ############################################
    pie_df = {'dummy':[], 'annotation':[], 'count':[]}
    for ann in annotation_classes:
        pie_df['dummy'].append('.')
        pie_df['annotation'].append(annotation_labels[ann])
        pie_df['count'].append(annotation_features[ann])

    ggplot.plot('%s/annotation_pie_pie.r'%os.environ['RDIR'], pie_df, [options.title, '%s_pie.pdf'%options.output_prefix])

    ############################################
    # read:length ratio
    ############################################
    ratio_df = {'annotation':[], 'ratio':[]}
    for ann in annotation_classes:
        ratio_df['annotation'].append(annotation_labels[ann])
        ratio_df['ratio'].append(annotation_ratio[ann])

    ggplot.plot('%s/annotation_pie_ratios.r'%os.environ['RDIR'], ratio_df, [options.title, '%s_ratios.pdf'%options.output_prefix])


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
    count = 0
    p = subprocess.Popen('intersectBed -s -u -f 0.5 -a %s -b %s' % (gff_file,bed_file), shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        count += 1
    p.communicate()
    return count


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
