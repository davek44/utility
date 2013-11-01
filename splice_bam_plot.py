#!/usr/bin/env python
from optparse import OptionParser
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.lib.ggplot2 as ggplot2
import math, pdb, subprocess
import gff

grdevices = importr('grDevices')

################################################################################
# splice_bam_plot.py
#
# Plot read coverage in a BAM file surrounding the splice sites defined in a
# gtf file.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf file> <bam file>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='control_bam_file', default=None, help='Control BAM file')
    parser.add_option('-o', dest='out_prefix', default='splice', help='Output prefix [Default: %default]')
    parser.add_option('-w', dest='window', default=200, type='int', help='Size of the surrounding window to consider [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide gtf file and BAM file')
    else:
        gtf_file = args[0]
        bam_file = args[1]

    # collect intervals
    intervals_5p, intervals_3p = get_splice_intervals(gtf_file, options.window)

    # process bam
    cov_5p = process_bam(bam_file, intervals_5p, options.out_prefix, options.window)
    cov_3p = process_bam(bam_file, intervals_3p, options.out_prefix, options.window)

    if options.control_bam_file:
        control_cov_5p = process_bam(options.control_bam_file, intervals_5p, options.out_prefix, options.window)
        control_cov_3p = process_bam(options.control_bam_file, intervals_3p, options.out_prefix, options.window)
    
    ############################################
    # output
    ############################################
    make_output(cov_5p, options.out_prefix+'_5p', options.window)
    make_output(cov_3p, options.out_prefix+'_3p', options.window)

    if options.control_bam_file:
        # normalize
        main_aligns = float(subprocess.check_output('samtools view -c %s' % bam_file, shell=True))
        control_aligns = float(subprocess.check_output('samtools view -c %s' % options.control_bam_file, shell=True))

        count_5p = sum([len(intervals_5p[chrom]) for chrom in intervals_5p])
        count_3p = sum([len(intervals_3p[chrom]) for chrom in intervals_3p])

        cov_5p_norm = [1000000.0*c/count_5p/main_aligns for c in cov_5p]
        cov_3p_norm = [1000000.0*c/count_3p/main_aligns for c in cov_3p]

        control_cov_5p_norm = [1000000.0*c/count_5p/control_aligns for c in control_cov_5p]
        control_cov_3p_norm = [1000000.0*c/count_3p/control_aligns for c in control_cov_3p]

        # plot and
        make_output_and(cov_5p_norm, control_cov_5p_norm, options.out_prefix+'_5p_and', options.window)
        make_output_and(cov_3p_norm, control_cov_3p_norm, options.out_prefix+'_3p_and', options.window)

        # plot subtraction
        sub_cov_5p = [cov_5p_norm[i]-control_cov_5p_norm[i] for i in range(len(cov_5p_norm))]
        sub_cov_3p = [cov_3p_norm[i]-control_cov_3p_norm[i] for i in range(len(cov_3p_norm))]

        make_output(sub_cov_5p, options.out_prefix+'_5p_sub', options.window)
        make_output(sub_cov_3p, options.out_prefix+'_3p_sub', options.window)


################################################################################
# get_splice_intervals
#
# Make a hash keyed by chromosome and valued with lists of (start,end,strand)
# TSS tuples. Sort by chromosome.
################################################################################
def get_splice_intervals(gtf_file, window):
    intervals_5p = {}
    intervals_3p = {}

    genes = gff.read_genes(gtf_file)
    for gid in genes:
        g = genes[gid]

        if len(g.exons) > 1:
            if g.strand == '+':
                # add first 5p site
                exon = g.exons[0]
                if not g.chrom in intervals_5p:
                    intervals_5p[g.chrom] = set()
                intervals_5p[g.chrom].add((exon.end - window/2, exon.end + window/2, g.strand))

                # process internal exons
                for exon in g.exons[1:-1]:
                    # add 3p site
                    if not g.chrom in intervals_3p:
                        intervals_3p[g.chrom] = set()
                    intervals_3p[g.chrom].add((exon.start - window/2, exon.start + window/2, g.strand))

                    # add 5p site
                    if not g.chrom in intervals_5p:
                        intervals_5p[g.chrom] = set()
                    intervals_5p[g.chrom].add((exon.end - window/2, exon.end + window/2, g.strand))

                # add last 3p site
                exon = g.exons[-1]
                if not g.chrom in intervals_3p:
                    intervals_3p[g.chrom] = set()
                intervals_3p[g.chrom].add((exon.start - window/2, exon.start + window/2, g.strand))

            else:
                # add first 5p site
                exon = g.exons[-1]
                if not g.chrom in intervals_5p:
                    intervals_5p[g.chrom] = set()
                intervals_5p[g.chrom].add((exon.end - window/2, exon.end + window/2, g.strand))

                # process internal exons (in reverse order, but doesn't matter)
                for exon in g.exons[1:-1]:
                    # add 3p site
                    if not g.chrom in intervals_3p:
                        intervals_3p[g.chrom] = set()
                    intervals_3p[g.chrom].add((exon.start - window/2, exon.start + window/2, g.strand))

                    # add 5p site
                    if not g.chrom in intervals_5p:
                        intervals_5p[g.chrom] = set()
                    intervals_5p[g.chrom].add((exon.end - window/2, exon.end + window/2, g.strand))

                # add last 3p site
                exon = g.exons[0]
                if not g.chrom in intervals_3p:
                    intervals_3p[g.chrom] = set()
                intervals_3p[g.chrom].add((exon.start - window/2, exon.start + window/2, g.strand))

    # convert sets to sorted lists
    for chrom in intervals_5p:
        intervals_5p[chrom] = sorted(list(intervals_5p[chrom]))
        intervals_3p[chrom] = sorted(list(intervals_3p[chrom]))

    return intervals_5p, intervals_3p


################################################################################
# make_bed
#
# Make a bed file of the TSS intervals to give to samtools
################################################################################
def make_bed(tss_intervals, bed_file):
    bed_out = open(bed_file,'w')
    for chrom in tss_intervals:
        for (start,end,strand) in tss_intervals[chrom]:
            print >> bed_out, '%s\t%d\t%d' % (chrom,start-1,end)
    bed_out.close()


################################################################################
# make_output
################################################################################
def make_output(cov, out_prefix, window):
    # dump raw counts to file
    raw_out = open('%s_raw.txt' % out_prefix,'w')
    for i in range(-window/2,window/2+1):
        print >> raw_out, '%d\t%e' % (i, cov[window/2+i])
    raw_out.close()

    # make plot data structures
    splice_i = ro.IntVector(range(-window/2,window/2+1))
    cov = ro.FloatVector(cov)
    df = ro.DataFrame({'splice_i':splice_i, 'cov':cov})

    # construct plot
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='splice_i', y='cov') + \
        ggplot2.geom_point() + \
        ggplot2.scale_x_continuous('Position relative to splice site') + \
        ggplot2.scale_y_continuous('Coverage')

    # plot to file
    grdevices.pdf(file='%s.pdf' % out_prefix)
    gp.plot()
    grdevices.dev_off()


################################################################################
# make_output_and
################################################################################
def make_output_and(cov, control_cov, out_prefix, window):
    # dump raw counts to file
    raw_out = open('%s_raw.txt' % out_prefix,'w')
    for i in range(-window/2,window/2+1):
        print >> raw_out, '%d\t%e\t%e' % (i, cov[window/2+i], control_cov[window/2+i])
    raw_out.close()

    # make plot data structures
    splice_i = ro.IntVector(2*range(-window/2,window/2+1))
    cov_r = ro.FloatVector(cov+control_cov)
    labels = ro.StrVector(['Main']*len(cov)+['Control']*len(control_cov))
    df = ro.DataFrame({'splice_i':splice_i, 'cov':cov_r, 'label':labels})

    # construct plot
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='splice_i', y='cov', colour='label') + \
        ggplot2.geom_point() + \
        ggplot2.scale_x_continuous('Position relative to splice site') + \
        ggplot2.scale_y_continuous('Coverage') + \
        ggplot2.scale_colour_discrete('')

    # plot to file
    grdevices.pdf(file='%s.pdf' % out_prefix)
    gp.plot()
    grdevices.dev_off()


################################################################################
# process_bam
#
# Count read coverage in a BAM file around the intervals given.
################################################################################
def process_bam(bam_file, intervals, out_prefix, window):
    # initialize data structures
    ti_chrom = None
    ti_i = -1
    active_intervals = []
    interval_cov = [0]*(window+1)

    make_bed(intervals,'%s_tss.bed' % out_prefix)

    proc = subprocess.Popen('samtools mpileup -l %s_tss.bed %s' % (out_prefix,bam_file), shell=True, stdout=subprocess.PIPE)
    line = proc.stdout.readline()
    while line:
        a = line.split()

        chrom = a[0]
        pos = int(a[1])
        cov = int(a[3])

        if chrom in intervals:
            # update active TSS
            active_intervals, ti_chrom, ti_i = update_active_intervals(intervals, ti_chrom, ti_i, active_intervals, chrom, pos)

            # append to interval_cov
            for tss in active_intervals:
                (tchrom,tstart,tend,tstrand) = tss
                if tstrand == '+':
                    interval_cov[pos-tstart] += cov
                else:
                    interval_cov[tend-pos] += cov

        # next line
        line = proc.stdout.readline()

    # finish
    line = proc.communicate()[0]

    return interval_cov


################################################################################
# update_active_intervals
#
# Maintain a set of TSS intervals which chrom and pos fall inside.
################################################################################
def update_active_intervals(intervals, ti_chrom, ti_i, active_intervals, chrom, pos):
    # drop finished tss
    while active_intervals and (active_intervals[0][0] != chrom or active_intervals[0][2] < pos):
        active_intervals = active_intervals[1:]

    # initialize intervals indexes
    if ti_chrom != chrom:
        ti_chrom = chrom
        ti_i = 0

    # add opened tss
    while ti_i < len(intervals[ti_chrom]):
        (tstart,tend,tstrand) = intervals[ti_chrom][ti_i]
        if pos < tstart:
            # before TSS: no new active
            break
        elif tstart <= pos <= tend:
            # inside TSS: add and move to next
            active_intervals.append((ti_chrom,tstart,tend,tstrand))
            ti_i += 1
        else:
            # past TSS: move to next
            ti_i += 1

    return active_intervals, ti_chrom, ti_i


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
