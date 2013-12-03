#!/usr/bin/env python
from optparse import OptionParser
import math, os, pdb, random, shutil, stats, subprocess, sys, tempfile
import pysam
import count_reads, gff, ggplot

################################################################################
# bam_heat.py
#
# Plot read coverage in a BAM file surrounding the median points of GFF entries
# as a heatmap.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gff> <bam>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='control_bam_file', default=None, help='Control BAM file')
    parser.add_option('-l', dest='log', default=False, action='store_true', help='log2 coverage [Default: %default]')
    parser.add_option('-m', dest='max_features', default=2000, type='int', help='Maximum number of features to plot [Default: %default]')
    parser.add_option('-o', dest='output_pre', default='bam', help='Output prefix [Default: %default]')
    parser.add_option('-s', dest='sorted', help='Plot heatmap in order of GFF [Default: %default]')
    parser.add_option('-u', dest='range', default=600, type='int', help='Range around peak middle [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide gtf file and BAM file')
    else:
        gff_file = args[0]
        bam_file = args[1]

    ############################################
    # extend GFF entries to range (and sample)
    ############################################
    feature_count = 0
    for line in open(gff_file):
        feature_count += 1

    sample_prob = min(1.0, options.max_features / float(feature_count))

    gff_range_fd, gff_range_file = tempfile.mkstemp()
    gff_range_out = open(gff_range_file, 'w')

    for line in open(gff_file):
        a = line.split('\t')
        
        start = int(a[3])
        end = int(a[4])
        mid = start + (end-start)/2
        a[3] = str(mid - options.range/2)
        a[4] = str(mid + options.range/2)

        if random.random() < sample_prob:
            print >> gff_range_out, '\t'.join(a),

    gff_range_out.close()

    ############################################
    # compute coverage
    ############################################
    coverage = compute_coverage(gff_range_file, bam_file)
    if options.control_bam_file:
        coverage_control = compute_coverage(gff_range_file, options.control_bam_file)

    ############################################
    # normalize
    ############################################
    # compute total fragments
    fragments = float(count_reads.count(bam_file, filter_mapq=True))
    if options.control_bam_file:
        fragments_control = float(count_reads.count(options.control_bam_file, filter_mapq=True))

    # normalize coverages
    for gff_entry in coverage:
        for i in range(len(coverage[gff_entry])):
            coverage[gff_entry][i] /= fragments
            if options.control_bam_file:
                coverage_control[gff_entry][i] /= fragments_control

    ############################################
    # plot heatmap
    ############################################
    if options.sorted:
        # gff file was sorted
        gff_entries_sorted = []
        for line in open(gff_range_file):
            a = line.split('\t')
            gchrom = a[0]
            gstart = int(a[3])
            gend = int(a[4])
            gff_entries_sorted.append((gchrom,gstart,gend))

    else:
        # sort by mean
        gff_entry_stat = []
        for gff_entry in coverage:
            if options.control_bam_file:
                gff_stat = stats.mean([math.log(1+coverage[gff_entry][i],2) - math.log(1+coverage_control[gff_entry][i],2) for i in range(len(coverage[gff_entry]))])
            else:
                gff_stat = stats.mean([coverage[gff_entry][i] for i in range(len(coverage[gff_entry]))])

            gff_entry_stat.append((gff_stat,gff_entry))

        gff_entry_stat.sort(reverse=True)

        gff_entries_sorted = [gff_entry for (gff_stat, gff_entry) in gff_entry_stat]

    df = {'Index':[], 'Feature':[], 'Coverage':[]}
    for g in range(len(gff_entries_sorted)):
        gff_entry = gff_entries_sorted[g]
        for i in range(-options.range/2,options.range/2+1):
            df['Index'].append(i)
            df['Feature'].append(g)

            if options.log:
                cov = math.log(1+coverage[gff_entry][-options.range/2+i],2)
            else:
                cov = coverage[gff_entry][-options.range/2+i]

            if options.control_bam_file:
                if options.log:
                    cov -= math.log(1+coverage_control[gff_entry][-options.range/2+i],2)
                else:
                    cov = (1+cov) / (1+coverage_control[gff_entry][-options.range/2+i])

            df['Coverage'].append('%.4e' % cov)

    r_script = '%s/bam_heat_heat.r' % os.environ['RDIR']
    out_pdf = '%s_heat.pdf' % options.output_pre

    ggplot.plot(r_script, df, [out_pdf], debug=True)

    ############################################
    # plot meta-coverage
    ############################################
    df = {'Index':[], 'Coverage':[]}
    if options.control_bam_file:
        df['Type'] = []

    for i in range(-options.range/2,options.range/2+1):
        df['Index'].append(i)

        if options.log:
            df['Coverage'].append(stats.geo_mean([(1+coverage[gff_entry][-options.range/2+i]) for gff_entry in coverage]))
        else:
            df['Coverage'].append(stats.mean([coverage[gff_entry][-options.range/2+i] for gff_entry in coverage]))

        if options.control_bam_file:
            df['Type'].append('Primary')

            df['Index'].append(i)
            df['Type'].append('Control')
            if options.log:
                df['Coverage'].append(stats.geo_mean([(1+coverage_control[gff_entry][-options.range/2+i]) for gff_entry in coverage_control]))
            else:
                df['Coverage'].append(stats.mean([coverage_control[gff_entry][-options.range/2+i] for gff_entry in coverage_control]))

    r_script = '%s/bam_heat_meta.r' % os.environ['RDIR']
    out_pdf = '%s_meta.pdf' % options.output_pre

    ggplot.plot(r_script, df, [out_pdf])


    # clean
    os.close(gff_range_fd)
    os.remove(gff_range_file)


################################################################################
# compute_coverage
#
# Input:
#  gff_file: GFF file of equal-sized genome features.
#  bam_file: BAM file of reads alignments.
################################################################################
def compute_coverage(gff_file, bam_file):
    # filter BAM for mapping quality
    bam_mapq_fd, bam_mapq_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
    bam_in = pysam.Samfile(bam_file, 'rb')
    bam_mapq_out = pysam.Samfile(bam_mapq_file, 'wb', template=bam_in)
    for aligned_read in bam_in:
        if aligned_read.mapq > 0:
            bam_mapq_out.write(aligned_read)
    bam_mapq_out.close()

    # count fragments and hash multi-mappers
    num_fragments = 0
    multi_maps = {}
    for aligned_read in pysam.Samfile(bam_mapq_file, 'rb'):
        try:
            nh_tag = aligned_read.opt('NH')
        except:
            nh_tag = 1.0

        if aligned_read.is_paired:
            num_fragments += 0.5/nh_tag
        else:
            num_fragments += 1.0/nh_tag

        if nh_tag > 1:
            multi_maps[aligned_read.qname] = nh_Tag

    # initialize counters
    coverage = {}
    for line in open(gff_file):
        a = line.split('\t')

        gchrom = a[0]
        gstart = int(a[3])
        gend = int(a[4])
        
        instance_id = (gchrom,gstart,gend)
        coverage[instance_id] = [0]*(gend-gstart+1)

    # count reads
    p = subprocess.Popen('intersectBed -split -wo -bed -abam %s -b %s' % (bam_mapq_file, gff_file), shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        a = line.split('\t')
        
        rstart = int(a[1])
        rend = int(a[2])
        rheader = a[3]

        # because intersectBed screws up indels near endpoints
        if rstart < rend:
            gchrom = a[12]
            gstart = int(a[15])
            gend = int(a[16])

            instance_id = (gchrom,gstart,gend)

            cov_start = max(rstart, gstart)
            cov_end = min(rend, gend)

            for i in range(cov_start - gstart, cov_end - gstart + 1):
                coverage[instance_id][i] += 1.0/multi_maps.get(rheader,1)

    p.communicate()

    # clean
    os.close(bam_mapq_fd)
    os.remove(bam_mapq_file)

    return coverage


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
