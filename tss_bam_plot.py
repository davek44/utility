#!/usr/bin/env python
from optparse import OptionParser
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.lib.ggplot2 as ggplot2
import math, pdb, subprocess
import gff

grdevices = importr('grDevices')

################################################################################
# tss_bam_plot.py
#
# Plot read coverage in a BAM file surrounding the TSS's defined in a gtf file.
#
# The geometric mean idea is nice in theory to smooth outliers, but the IP
# sequencing will have more variance than the Control which throws off the
# the alignment # normalization.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf file> <bam file>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='control_bam_file', default=None, help='Control BAM file')
    parser.add_option('-d', dest='downstream', default=2000, type='int', help='TSS downstream [Default: %default]')
    parser.add_option('-g', dest='geo_mean', default=False, action='store_true', help='Plot coverage geometric means [Default: %default]')
    parser.add_option('-o', dest='out_prefix', default='tss', help='Output prefix [Default: %default]')
    parser.add_option('-u', dest='upstream', default=5000, type='int', help='TSS upstream [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide gtf file and BAM file')
    else:
        gtf_file = args[0]
        bam_file = args[1]

    # collect intervals
    tss_intervals, tss_count = get_tss(gtf_file, options.upstream, options.downstream)

    # process bam
    tss_cov = process_bam(bam_file, tss_intervals, options.out_prefix, options.upstream, options.downstream, options.geo_mean)

    if options.control_bam_file:
        control_tss_cov = process_bam(options.control_bam_file, tss_intervals, options.out_prefix, options.upstream, options.downstream, options.geo_mean)
    
    ############################################
    # output
    ############################################
    make_output(tss_cov, options.out_prefix, options.upstream, options.downstream)

    if options.control_bam_file:
        # normalize
        main_aligns = float(subprocess.check_output('samtools view -c %s' % bam_file, shell=True))
        control_aligns = float(subprocess.check_output('samtools view -c %s' % options.control_bam_file, shell=True))

        if options.geo_mean:
            tss_cov = [1000000.0*math.exp(float(tc)/tss_count)/main_aligns for tc in tss_cov]
            control_tss_cov = [1000000.0*math.exp(float(tc)/tss_count)/control_aligns for tc in control_tss_cov]
        else:
            tss_cov = [1000000.0*tc/tss_count/main_aligns for tc in tss_cov]
            control_tss_cov = [1000000.0*tc/tss_count/control_aligns for tc in control_tss_cov]

        # plot subtraction
        sub_cov = [tss_cov[i]-control_tss_cov[i] for i in range(len(tss_cov))]
        make_output(sub_cov, options.out_prefix+'_sub', options.upstream, options.downstream)

        # plot and
        make_output_and(tss_cov, control_tss_cov, options.out_prefix+'_and', options.upstream, options.downstream)


################################################################################
# get_tss
#
# Make a hash keyed by chromosome and valued with lists of (start,end,strand)
# TSS tuples. Sort by chromosome.
################################################################################
def get_tss(gtf_file, upstream, downstream):
    tss_intervals = {}
    tss_count = 0

    genes = gff.read_genes(gtf_file)
    for gid in genes:
        g = genes[gid]
        if g.strand == '+':
            istart = g.exons[0].start-upstream
            iend = g.exons[0].start+downstream
        else:
            istart = g.exons[-1].end-downstream
            iend = g.exons[-1].end+upstream

        tss_intervals.setdefault(g.chrom,[]).append((istart,iend,g.strand))
        tss_count += 1

    for chrom in tss_intervals:
        tss_intervals[chrom].sort()

    return tss_intervals, tss_count


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
def make_output(tss_cov, out_prefix, upstream, downstream):
    # dump raw counts to file
    raw_out = open('%s_raw.txt' % out_prefix,'w')
    for i in range(-upstream,downstream+1):
        print >> raw_out, '%d\t%e' % (i, tss_cov[upstream+i])
    raw_out.close()

    # make plot data structures
    tss_i = ro.IntVector(range(-upstream,downstream+1))
    cov = ro.FloatVector(tss_cov)
    df = ro.DataFrame({'tss_i':tss_i, 'cov':cov})

    # construct full plot
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='tss_i', y='cov') + \
        ggplot2.geom_point() + \
        ggplot2.scale_x_continuous('TSS index') + \
        ggplot2.scale_y_continuous('Coverage')

    # plot to file
    grdevices.pdf(file='%s_full.pdf' % out_prefix)
    gp.plot()
    grdevices.dev_off()

    # construct zoomed plot
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='tss_i', y='cov') + \
        ggplot2.geom_point() + \
        ggplot2.scale_x_continuous('TSS index',limits=ro.IntVector([-1000,1000])) + \
        ggplot2.scale_y_continuous('Coverage')

    # plot to file
    grdevices.pdf(file='%s_zoom.pdf' % out_prefix)
    gp.plot()
    grdevices.dev_off()


################################################################################
# make_output_and
################################################################################
def make_output_and(tss_cov, control_tss_cov, out_prefix, upstream, downstream):
    # dump raw counts to file
    raw_out = open('%s_raw.txt' % out_prefix,'w')
    for i in range(-upstream,downstream+1):
        print >> raw_out, '%d\t%e\t%e' % (i, tss_cov[upstream+i], control_tss_cov[upstream+i])
    raw_out.close()

    # make plot data structures
    tss_i = ro.IntVector(2*range(-upstream,downstream+1))
    cov = ro.FloatVector(tss_cov+control_tss_cov)
    labels = ro.StrVector(['Main']*len(tss_cov)+['Control']*len(control_tss_cov))
    df = ro.DataFrame({'tss_i':tss_i, 'cov':cov, 'label':labels})

    # construct full plot
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='tss_i', y='cov', colour='label') + \
        ggplot2.geom_point() + \
        ggplot2.scale_x_continuous('TSS index') + \
        ggplot2.scale_y_continuous('Coverage') + \
        ggplot2.scale_colour_discrete('')

    # plot to file
    grdevices.pdf(file='%s_full.pdf' % out_prefix)
    gp.plot()
    grdevices.dev_off()

    # construct zoomed plot
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='tss_i', y='cov', colour='label') + \
        ggplot2.geom_point() + \
        ggplot2.scale_x_continuous('TSS index',limits=ro.IntVector([-1000,1000])) + \
        ggplot2.scale_y_continuous('Coverage') + \
        ggplot2.scale_colour_discrete('')

    # plot to file
    grdevices.pdf(file='%s_zoom.pdf' % out_prefix)
    gp.plot()
    grdevices.dev_off()


################################################################################
# process_bam
#
# Count read coverage in a BAM file around the tss_intervals given.
################################################################################
def process_bam(bam_file, tss_intervals, out_prefix, upstream, downstream, geo_mean):
    # initialize data structures
    ti_chrom = None
    ti_i = -1
    active_tss = []
    tss_cov = [0]*(upstream+downstream+1)

    make_bed(tss_intervals,'%s_tss.bed' % out_prefix)

    proc = subprocess.Popen('samtools mpileup -l %s_tss.bed %s' % (out_prefix,bam_file), shell=True, stdout=subprocess.PIPE)
    line = proc.stdout.readline()
    while line:
        a = line.split()

        chrom = a[0]
        pos = int(a[1])
        cov = int(a[3])

        if chrom in tss_intervals:
            # update active TSS
            active_tss, ti_chrom, ti_i = update_active_tss(tss_intervals, ti_chrom, ti_i, active_tss, chrom, pos)

            # append to tss_cov
            for tss in active_tss:
                (tchrom,tstart,tend,tstrand) = tss
                if tstrand == '+':
                    if geo_mean:
                        tss_cov[pos-tstart] += math.log(cov)
                    else:
                        tss_cov[pos-tstart] += cov
                else:
                    if geo_mean:
                        tss_cov[tend-pos] += math.log(cov)
                    else:
                        tss_cov[tend-pos] += cov

        # next line
        line = proc.stdout.readline()

    # finish
    line = proc.communicate()[0]

    return tss_cov


################################################################################
# update_active_tss
#
# Maintain a set of TSS intervals which chrom and pos fall inside.
################################################################################
def update_active_tss(tss_intervals, ti_chrom, ti_i, active_tss, chrom, pos):
    # drop finished tss
    while active_tss and (active_tss[0][0] != chrom or active_tss[0][2] < pos):
        active_tss = active_tss[1:]

    # initialize tss_intervals indexes
    if ti_chrom != chrom:
        ti_chrom = chrom
        ti_i = 0

    # add opened tss
    while ti_i < len(tss_intervals[ti_chrom]):
        (tstart,tend,tstrand) = tss_intervals[ti_chrom][ti_i]
        if pos < tstart:
            # before TSS: no new active
            break
        elif tstart <= pos <= tend:
            # inside TSS: add and move to next
            active_tss.append((ti_chrom,tstart,tend,tstrand))
            ti_i += 1
        else:
            # past TSS: move to next
            ti_i += 1

    return active_tss, ti_chrom, ti_i


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
