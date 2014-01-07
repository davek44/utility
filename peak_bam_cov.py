#!/usr/bin/env python
from optparse import OptionParser
import os, pdb, shutil, stats, subprocess, sys, tempfile
import pysam
import gff, ggplot

################################################################################
# peak_bam_cov.py
#
# Plot read coverage in a BAM file surrounding the median points of GFF entries
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gff> <bam>'
    parser = OptionParser(usage)
    #parser.add_option('-c', dest='control_bam_file', default=None, help='Control BAM file')
    parser.add_option('-g', dest='geo_mean', default=False, action='store_true', help='Compute geometric mean of individual peak coverages [Default: %default]')
    parser.add_option('-i', dest='individual_plots', default=False, action='store_true', help='Print a coverage plot for every individual peak [Default: %default]')
    parser.add_option('-o', dest='out_prefix', default='peak_cov', help='Output prefix [Default: %default]')
    parser.add_option('-p', dest='properly_paired', default=False, action='store_true', help='Count entire fragments for only properly paired reads [Default: %default]')
    parser.add_option('-u', dest='range', default=500, type='int', help='Range around peak middle [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide gtf file and BAM file')
    else:
        peaks_gff = args[0]
        bam_file = args[1]

    # filter BAM for mapping quality
    bam_mapq_fd, bam_mapq_file = tempfile.mkstemp(dir='%s/research/scratch' % os.environ['HOME'])
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
        if options.properly_paired:
            if aligned_read.is_properly_paired:
                num_fragments += 0.5/aligned_read.opt('NH')
        else:
            if aligned_read.is_paired:
                num_fragments += 0.5/aligned_read.opt('NH')
            else:
                num_fragments += 1.0/aligned_read.opt('NH')

        if aligned_read.opt('NH') > 1:
            multi_maps[aligned_read.qname] = aligned_read.opt('NH')

    # extend GFF entries to range
    peaks_gff_range_fd, peaks_gff_range_file = tempfile.mkstemp()
    peaks_gff_range_out = open(peaks_gff_range_file, 'w')
    for line in open(peaks_gff):
        a = line.split('\t')
        
        pstart = int(a[3])
        pend = int(a[4])
        peak_mid = pstart + (pend-pstart)/2

        a[3] = str(peak_mid - options.range/2 - 1)
        a[4] = str(peak_mid + options.range/2 + 1)

        print >> peaks_gff_range_out, '\t'.join(a),
    peaks_gff_range_out.close()

    # initialize coverage counters
    peak_cov_individual = {}
    peak_reads = {}

    # count reads
    p = subprocess.Popen('intersectBed -split -wo -bed -abam %s -b %s' % (bam_mapq_file,peaks_gff_range_file), shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        a = line.split('\t')
        
        rstart = int(a[1])
        rend = int(a[2])
        rheader = a[3]

        # because intersectBed screws up indels near endpoints
        if rstart < rend:
            pstart = int(a[15])
            pend = int(a[16])
            peak_id = gff.gtf_kv(a[20])['id']
            peak_reads[peak_id] = peak_reads.get(peak_id,0) + 1

            peak_mid = pstart + (pend-pstart)/2
            peak_range_start = peak_mid - options.range/2
            peak_range_end = peak_mid + options.range/2

            range_start = max(rstart, peak_range_start)
            range_end = min(rend, peak_range_end)

            if not peak_id in peak_cov_individual:
                peak_cov_individual[peak_id] = [0.0]*(1+options.range)
            for i in range(range_start - peak_range_start, range_end - peak_range_start + 1):
                peak_cov_individual[peak_id][i] += 1.0/multi_maps.get(rheader,1)

    p.communicate()

    # combine individual
    peak_cov = [0.0]*(1+options.range)
    for i in range(len(peak_cov)):
        if options.geo_mean:
            peak_cov[i] = stats.geo_mean([1+peak_cov_individual[peak_id][i] for peak_id in peak_cov_individual])
        else:
            peak_cov[i] = stats.mean([peak_cov_individual[peak_id][i] for peak_id in peak_cov_individual])

    #for peak_id in peak_reads:
    #    print peak_id, peak_reads[peak_id]

    # output
    make_output(peak_cov, options.out_prefix, options.range)

    if options.individual_plots:
        individual_dir = '%s_individuals' % options.out_prefix
        if os.path.isdir(individual_dir):
            shutil.rmtree(individual_dir)
        os.mkdir(individual_dir)

        for peak_id in peak_cov_individual:
            if peak_reads[peak_id] > 150:
                make_output(peak_cov_individual[peak_id], '%s/%s' % (individual_dir,peak_id), options.range)

    # clean
    os.close(bam_mapq_fd)
    os.remove(bam_mapq_file)
    os.close(peaks_gff_range_fd)
    os.remove(peaks_gff_range_file)


################################################################################
# make_output
################################################################################
def make_output(peak_cov, out_prefix, prange):
    # dump raw counts to file
    raw_out = open('%s_raw.txt' % out_prefix,'w')
    for i in range(-prange/2,prange/2+1):
        print >> raw_out, '%d\t%e' % (i, peak_cov[i+prange/2])
    raw_out.close()

    r_script = '%s/peak_bam_plot.r' % os.environ['RDIR']

    df_dict = {'peak_i':range(-prange/2,prange/2+1),
               'cov':peak_cov}

    out_pdf = '%s.pdf' % out_prefix

    ggplot.plot(r_script, df_dict, [out_pdf])


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
