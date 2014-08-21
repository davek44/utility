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
    usage = 'usage: %prog [options] <gff> <bam1,bam2,...>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='control_bam_files', default=None, help='Control BAM files (comma separated)')
    parser.add_option('-l', dest='log', default=False, action='store_true', help='log2 coverage [Default: %default]')
    parser.add_option('-k', dest='gtf_key', default=None, help='GTF key to hash gff entries by')
    parser.add_option('-m', dest='max_features', default=2000, type='int', help='Maximum number of features to plot [Default: %default]')
    parser.add_option('-o', dest='output_pre', default='bam', help='Output prefix [Default: %default]')
    parser.add_option('-s', dest='sorted_gene_files', help='Files of sorted gene lists. Plot heatmaps in their order')
    parser.add_option('-u', dest='range', default=2000, type='int', help='Range around peak middle [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide gtf file and BAM file')
    else:
        gff_file = args[0]
        bam_files = args[1].split(',')

    if options.control_bam_files:
        control_bam_files = options.control_bam_files.split(',')

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

        range_start = mid - options.range/2
        range_end = mid + options.range/2

        if range_start > 0:
            a[3] = str(mid - options.range/2)
            a[4] = str(mid + options.range/2)
            a[-1] = a[-1].rstrip()

            if random.random() < sample_prob:
                print >> gff_range_out, '\t'.join(a)

    gff_range_out.close()

    ############################################
    # compute coverage
    ############################################
    coverage, fragments = compute_coverage(gff_range_file, bam_files, options.gtf_key)
    if options.control_bam_files:
        coverage_control, fragments_control = compute_coverage(gff_range_file, control_bam_files, options.gtf_key)

    # clean
    os.close(gff_range_fd)
    os.remove(gff_range_file)

    ############################################
    # normalize
    ############################################
    # normalize coverages (and add pseudocounts)
    for feature_id in coverage:
        for i in range(len(coverage[feature_id])):
            coverage[feature_id][i] = (1+coverage[feature_id][i])/fragments
            if options.control_bam_files:
                coverage_control[feature_id][i] = (1+coverage_control[feature_id][i])/fragments_control    

    ############################################
    # sorted genes
    ############################################
    features_sorted = []
    if options.sorted_gene_files:
        # for each sorted list
        for sorted_gene_file in options.sorted_gene_files.split(','):
            # collect feature_id's
            features_sorted.append([])
            for line in open(sorted_gene_file):
                feature_id = line.split()[0]
                # verify randomly selected
                if feature_id in coverage:
                    features_sorted[-1].append(feature_id)

    else:
        # tuple feature_id's with mean coverage
        feature_id_stat = []
        for feature_id in coverage:
            if options.control_bam_files:
                feature_stat = stats.mean([math.log(coverage[feature_id][i],2) - math.log(coverage_control[feature_id][i],2) for i in range(len(coverage[feature_id]))])
            else:
                feature_stat = stats.geo_mean([coverage[feature_id][i] for i in range(len(coverage[feature_id]))])

            feature_id_stat.append((feature_stat,feature_id))

        # sort
        feature_id_stat.sort(reverse=True)

        # store as the only sorted list
        features_sorted.append([feature_id for (feature_stat, feature_id) in feature_id_stat])

    ############################################
    # plot heatmap(s)
    ############################################
    # if multiple sorts, create a dir for the plots
    if len(features_sorted) > 1:
        if not os.path.isdir('%s_heat' % options.output_pre):
            os.mkdir('%s_heat' % options.output_pre)

    for s in range(len(features_sorted)):
        df = {'Index':[], 'Feature':[], 'Coverage':[]}
        for f in range(len(features_sorted[s])):
            feature_id = features_sorted[s][f]
            for i in range(-options.range/2,options.range/2+1):
                df['Index'].append(i)
                df['Feature'].append(f)

                if options.log:
                    cov = math.log(coverage[feature_id][i+options.range/2],2)
                else:
                    cov = coverage[feature_id][i+options.range/2]

                if options.control_bam_files:
                    if options.log:
                        cov -= math.log(coverage_control[feature_id][i+options.range/2],2)
                    else:
                        cov = cov / coverage_control[feature_id][i+options.range/2]

                df['Coverage'].append('%.4e' % cov)

        r_script = '%s/bam_heat_heat.r' % os.environ['RDIR']
        if len(features_sorted) == 1:
            out_pdf = '%s_heat.pdf' % options.output_pre
        else:
            sorted_gene_file = options.sorted_gene_files.split(',')[s]
            sorted_gene_pre = os.path.splitext(os.path.split(sorted_gene_file)[-1])[0]
            out_pdf = '%s_heat/%s.pdf' % (options.output_pre,sorted_gene_pre)

        ggplot.plot(r_script, df, [out_pdf, options.control_bam_files!=None])

    ############################################
    # plot meta-coverage
    ############################################
    df = {'Index':[], 'Coverage':[]}
    if options.control_bam_files:
        df['Type'] = []

    for i in range(-options.range/2,options.range/2+1):
        df['Index'].append(i)

        if options.log:
            df['Coverage'].append(stats.geo_mean([coverage[feature_id][i+options.range/2] for feature_id in coverage]))
        else:
            df['Coverage'].append(stats.mean([coverage[feature_id][i+options.range/2] for feature_id in coverage]))

        if options.control_bam_files:
            df['Type'].append('Primary')

            df['Index'].append(i)
            df['Type'].append('Control')
            if options.log:
                df['Coverage'].append(stats.geo_mean([coverage_control[feature_id][i+options.range/2] for feature_id in coverage_control]))
            else:
                df['Coverage'].append(stats.mean([coverage_control[feature_id][i+options.range/2] for feature_id in coverage_control]))

    r_script = '%s/bam_heat_meta.r' % os.environ['RDIR']
    out_pdf = '%s_meta.pdf' % options.output_pre

    ggplot.plot(r_script, df, [out_pdf])


################################################################################
# compute_coverage
#
# Input:
#  gff_file: GFF file of equal-sized genome features.
#  bam_file: BAM file of reads alignments.
#  gtf_key:  GTF key by which is hash coverage arrays.
################################################################################
def compute_coverage(gff_file, bam_files, gtf_key):
    # initialize counters
    fragments = 0
    coverage = {}
    for line in open(gff_file):
        a = line.split('\t')

        gchrom = a[0]
        gstart = int(a[3])
        gend = int(a[4])

        if gtf_key == None:
            instance_id = (gchrom,gstart,gend)
        else:
            instance_id = gff.gtf_kv(a[8])[gtf_key]

        coverage[instance_id] = [0]*(gend-gstart+1)

    # process bam files
    for bam_file in bam_files:
        # filter BAM for mapping quality
        bam_mapq_fd, bam_mapq_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
        bam_in = pysam.Samfile(bam_file, 'rb')
        bam_mapq_out = pysam.Samfile(bam_mapq_file, 'wb', template=bam_in)
        for aligned_read in bam_in:
            if aligned_read.mapq > 0:
                bam_mapq_out.write(aligned_read)
        bam_mapq_out.close()

        # count fragments and hash multi-mappers
        multi_maps = {}
        paired_reads = False
        for aligned_read in pysam.Samfile(bam_mapq_file, 'rb'):
            try:
                nh_tag = aligned_read.opt('NH')
            except:
                nh_tag = 1.0

            if aligned_read.is_paired:
                paired_reads = True
                fragments += 0.5/nh_tag
            else:
                fragments += 1.0/nh_tag

            if nh_tag > 1:
                multi_maps[aligned_read.qname] = nh_tag

        # count reads
        p = subprocess.Popen('intersectBed -split -wo -bed -abam %s -b %s' % (bam_mapq_file, gff_file), shell=True, stdout=subprocess.PIPE)
        for line in p.stdout:
            a = line.split('\t')

            rstart = int(a[1])+1  # convert back to 1-based
            rend = int(a[2])
            rheader = a[3]

            # because intersectBed screws up indels near endpoints
            if rstart < rend:
                gchrom = a[12]
                gstart = int(a[15])
                gend = int(a[16])
                gstrand = a[18]

                if gtf_key == None:
                    instance_id = (gchrom,gstart,gend)
                else:
                    instance_id = gff.gtf_kv(a[20])[gtf_key]

                cov_start = max(rstart, gstart)
                cov_end = min(rend, gend)

                if gstrand == '+':
                    inc_start = cov_start - gstart
                    inc_end = cov_end - gstart + 1
                else:
                    inc_start = gend - cov_end
                    inc_end = gend - cov_start + 1

                # find multi-map number, which may require removing a suffix
                if rheader in multi_maps:
                    mm = multi_maps[rheader]
                else:
                    rheader_base = rheader[:rheader.rfind('/')]
                    if rheader_base in multi_maps:
                        mm = multi_maps[rheader_base]
                    else:
                        mm = 1.0

                for i in range(inc_start, inc_end):
                    coverage[instance_id][i] += 1.0/mm

        p.communicate()

        # clean
        os.close(bam_mapq_fd)
        os.remove(bam_mapq_file)

    return coverage, fragments


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
