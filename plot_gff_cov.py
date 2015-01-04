#!/usr/bin/env python
from optparse import OptionParser
import math, os, pdb, random, shutil, stats, subprocess, sys, tempfile

#from guppy import hpy
import pysam

import gff, ggplot

################################################################################
# plot_gff_cov.py
#
# Plot coverage of entries in a BAM or GFF file over the median points or span
# of GFF entries as a heatmap.
#
# To plot coverage around TSS:
#  ./plot_gff_cov.py -u mid -r 1000 mid tss.gff reads.bam
#
# To plot coverage across gene spans:
#  ./plot_gff_cov.py -m 500 span genes.gtf reads.bam
#
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <mode=mid/span> <anchor_gff> <event_bam1,event_bam2,...|event_gff1,event_gff2,...>'
    parser = OptionParser(usage)
    parser.add_option('-a', dest='max_anchors', default=1000, type='int', help='Maximum number of anchors to consider [Default: %default]')
    parser.add_option('-c', dest='control_files', default=None, help='Control BAM or GFF files (comma separated)')
    parser.add_option('-e', dest='plot_heat', default=False, action='store_true', help='Plot as a heatmap [Default: %default]')
    parser.add_option('-l', dest='log', default=False, action='store_true', help='log2 coverage [Default: %default]')
    parser.add_option('--labels', dest='labels', default='Primary,Control', help='Plot labels [Default:%default]')
    parser.add_option('-o', dest='output_pre', default='gff_cov', help='Output prefix [Default: %default]')
    parser.add_option('-s', dest='sorted_gene_files', help='Files of sorted gene lists. Plot heatmaps in their order')

    parser.add_option('-p', dest='smooth_span', default=0.2, type='float', help='Smoothing span parameter [Default: %default]')

    parser.add_option('-b', dest='bins', default=100, type='int', help='Number of bins across the gene span [Default: %default]')
    parser.add_option('-m', dest='min_length', default=None, type='int', help='Minimum anchor length [Default: %default]')

    parser.add_option('-w', dest='window', default=2000, type='int', help='Window around peak middle [Default: %default]')

    (options,args) = parser.parse_args()

    if len(args) != 3:
        parser.error('Must provide mode, anchor GFF, and BAM/GFF file(s)')
    else:
        mode = args[0]
        anchor_gff = args[1]
        event_files = args[2].split(',')

    if options.control_files:
        control_files = options.control_files.split(',')

    plot_labels = options.labels.split(',')

    anchor_is_gtf = (anchor_gff[-4:] == '.gtf')

    # preprocess anchor GFF
    prep_anchor_fd, prep_anchor_gff = preprocess_anchors(anchor_gff, mode, options.max_anchors, anchor_is_gtf, options.min_length, options.window)

    ############################################
    # compute coverage
    ############################################
    coverage, events = compute_coverage(prep_anchor_gff, event_files, mode, anchor_is_gtf, options.bins)
    if options.control_files:
        coverage_control, events_control = compute_coverage(prep_anchor_gff, control_files, mode, anchor_is_gtf, options.bins)

    # clean
    os.close(prep_anchor_fd)
    os.remove(prep_anchor_gff)

    ############################################
    # normalize
    ############################################
    # normalize coverages (and add pseudocounts)
    for anchor_id in coverage:
        for i in range(len(coverage[anchor_id])):
            coverage[anchor_id][i] = (1+coverage[anchor_id][i])/float(events)
            if options.control_files:
                coverage_control[anchor_id][i] = (1+coverage_control[anchor_id][i])/float(events_control)

    ############################################
    # sort anchors
    ############################################
    anchors_sorted = []
    if options.sorted_gene_files:
        # for each sorted list
        for sorted_gene_file in options.sorted_gene_files.split(','):
            # collect anchor_id's
            anchors_sorted.append([])
            for line in open(sorted_gene_file):
                anchor_id = line.split()[0]
                # verify randomly selected
                if anchor_id in coverage:
                    anchors_sorted[-1].append(anchor_id)

    else:
        # tuple anchor_id's with mean coverage
        stat_aid = []
        for anchor_id in coverage:
            if options.control_files:
                astat = stats.mean([math.log(coverage[anchor_id][i],2) - math.log(coverage_control[anchor_id][i],2) for i in range(len(coverage[anchor_id]))])
            else:
                astat = stats.geo_mean([coverage[anchor_id][i] for i in range(len(coverage[anchor_id]))])

            stat_aid.append((astat, anchor_id))

        # sort
        stat_aid.sort(reverse=True)

        # store as the only sorted list
        anchors_sorted.append([anchor_id for (astat, anchor_id) in stat_aid])

    ############################################
    # plot heatmap(s)
    ############################################
    if options.plot_heat:
        # if multiple sorts, create a dir for the plots
        if len(anchors_sorted) > 1:
            if not os.path.isdir('%s_heat' % options.output_pre):
                os.mkdir('%s_heat' % options.output_pre)

        for s in range(len(anchors_sorted)):
            df = {'Index':[], 'Anchor':[], 'Coverage':[]}
            for si in range(len(anchors_sorted[s])):
                anchor_id = anchors_sorted[s][si]

                for i in range(len(coverage[anchor_id])):
                    if mode == 'mid':
                        df['Index'].append(i - options.window/2)
                    else:
                        df['Index'].append(i)
                    df['Anchor'].append(anchor_id)

                    if options.log:
                        cov = math.log(coverage[anchor_id][i], 2)
                    else:
                        cov = coverage[anchor_id][i]

                    if options.control_files:
                        if options.log:
                            cov -= math.log(coverage_control[anchor_id][i], 2)
                        else:
                            cov = cov / coverage_control[anchor_id][i]

                    df['Coverage'].append('%.4e' % cov)
            
            if len(anchors_sorted) == 1:
                out_pdf = '%s_heat.pdf' % options.output_pre
            else:
                sorted_gene_file = options.sorted_gene_files.split(',')[s]
                sorted_gene_pre = os.path.splitext(os.path.split(sorted_gene_file)[-1])[0]
                out_pdf = '%s_heat/%s.pdf' % (options.output_pre,sorted_gene_pre)

            r_script = '%s/plot_gff_cov_heat.r' % os.environ['RDIR']
            ggplot.plot(r_script, df, [out_pdf, options.control_files!=None])

    ############################################
    # plot meta-coverage
    ############################################
    df = {'Index':[], 'Coverage':[]}
    if options.control_files:
        df['Type'] = []

    if mode == 'mid':
        index_length = 2*(options.window/2) + 1
    elif mode == 'span':
        index_length = options.bins
    else:
        print >> sys.stderr, 'Unknown mode %s' % mode
        exit(1)

    for i in range(index_length):
        if mode == 'mid':
            df['Index'].append(i - options.window/2)
        else:
            df['Index'].append(i)

        if options.log:
            df['Coverage'].append(stats.geo_mean([coverage[anchor_id][i] for anchor_id in coverage]))
        else:
            df['Coverage'].append(stats.mean([coverage[anchor_id][i] for anchor_id in coverage]))

        if options.control_files:
            df['Type'].append('Primary')

            if mode == 'mid':
                df['Index'].append(i - options.window/2)
            else:
                df['Index'].append(i)

            df['Type'].append('Control')
            if options.log:
                df['Coverage'].append(stats.geo_mean([coverage_control[anchor_id][i] for anchor_id in coverage_control]))
            else:
                df['Coverage'].append(stats.mean([coverage_control[anchor_id][i] for anchor_id in coverage_control]))

    r_script = '%s/plot_gff_cov_meta.r' % os.environ['RDIR']
    out_df = '%s_meta.df' % options.output_pre
    ggplot.plot(r_script, df, [options.output_pre, options.smooth_span, plot_labels[0], plot_labels[1]], df_file=out_df)


################################################################################
# compute_coverage
#
# Input
#  anchor_gff:    GFF file of equal-sized genome features.
#  event_files:   BAM or GFF files of reads alignments.
#  mode:          mid or span.
#  anchor_is_gtf: True iff anchor_gff is GTF.
#  bins:          Number of bins to consider in span mode.
#
# Output
#  coverage:      Dict mapping anchor_id's to coverage arrays.
#  events:        Total number of events.
################################################################################
def compute_coverage(anchor_gff, event_files, mode, anchor_is_gtf, bins):
    ############################################
    # initialize
    ############################################
    coverage = initialize_coverage(anchor_gff, mode, anchor_is_gtf, bins)    

    if anchor_is_gtf:
        # get transcript structures
        transcripts = gff.read_genes(anchor_gff, key_id='transcript_id')

        # compute lengths
        transcript_lengths = {}
        for tid in transcripts:
            tx = transcripts[tid]
            for exon in tx.exons:
                transcript_lengths[tid] = transcript_lengths.get(tid,0) + exon.end-exon.start+1

    else:
        transcripts = None
        transcript_lengths = None


    events = 0
    for event_file in event_files:
        print >> sys.stderr, 'Computing coverage for %s in %s' % (event_file, anchor_gff)

        ############################################
        # preprocess BAM/GFF
        ############################################
        if event_file[-4:] == '.bam':
            # count fragments and hash multi-mappers
            multi_maps = {}
            for aligned_read in pysam.Samfile(event_file, 'rb'):
                try:
                    nh_tag = aligned_read.opt('NH')
                except:
                    nh_tag = 1.0

                if aligned_read.is_paired:
                    events += 0.5/nh_tag
                else:
                    events += 1.0/nh_tag

                if nh_tag > 1:
                    multi_maps[aligned_read.qname] = nh_tag

        elif event_file[-4:] == '.gff':
            for line in open(event_file):
                events += 1

        else:
            print >> sys.stderr, 'Unknown event file format %s' % event_file

        ############################################
        # intersect BAM w/ anchors
        ############################################
        if event_file[-4:] == '.bam':
            p = subprocess.Popen('intersectBed -split -wo -bed -abam %s -b %s' % (event_file, anchor_gff), shell=True, stdout=subprocess.PIPE)
        else:
            p = subprocess.Popen('intersectBed -s -wo -a %s -b %s' % (event_file, anchor_gff), shell=True, stdout=subprocess.PIPE)

        for line in p.stdout:
            a = line.split('\t')

            if event_file[-4:] == '.bam':
                rstart = int(a[1])+1  # convert back to 1-based gff from bed
                rend = int(a[2])
                rheader = a[3]
            else:
                rstart = int(a[3])
                rend = int(a[4])

            # because intersectBed screws up indels near endpoints
            if rstart < rend:
                if event_file[-4:] == '.bam':
                    acol = 12
                else:
                    acol = 9

                achrom = a[acol]
                astart = int(a[acol+3])
                aend = int(a[acol+4])
                astrand = a[acol+6]

                if anchor_is_gtf:
                    anchor_id = gff.gtf_kv(a[acol+8])['transcript_id']
                else:
                    anchor_id = '%s:%d-%d' % (achrom, astart, aend)

                # find where to increment
                inc_start, inc_end = find_inc_coords(anchor_id, astart, aend, astrand, rstart, rend, mode, bins, transcripts, transcript_lengths)

                if inc_start != None:
                    if event_file[-4:] == '.bam':
                        # find multi-map number, which may require removing a suffix
                        if rheader in multi_maps:
                            mm = multi_maps[rheader]
                        else:
                            rheader_base = rheader[:rheader.rfind('/')]
                            if rheader_base in multi_maps:
                                mm = multi_maps[rheader_base]
                            else:
                                mm = 1.0
                    else:
                        mm = 1.0

                    # increment!
                    for i in range(inc_start, inc_end):
                        coverage[anchor_id][i] += 1.0/mm

        p.communicate()

    return coverage, events


################################################################################
# find_inc_coords
#
# Input
#  anchor_id:   Anchor ID.
#  astart:      Anchor start.
#  aend:        Anchor end.
#  astrand:     Anchor strand.
#  rstart:      Read alignment start.
#  rend:        Read alignment end.
#  mode:        mid or span.
#  bins:        Number of bins in span mode.
#  transcripts: Transcript structures.
#  transcript_lengths:
#
# Output
#  inc_start:   Coordinate at which to start incrementing.
#  inc_end:     Coordinate at which to stop incrementing.
################################################################################
def find_inc_coords(anchor_id, astart, aend, astrand, rstart, rend, mode, bins, transcripts, transcript_lengths):
    if mode == 'mid':
        cov_start = max(rstart, astart)
        cov_end = min(rend, aend)

        if astrand == '+':
            inc_start = cov_start - astart
            inc_end = cov_end - astart + 1
        else:
            inc_start = aend - cov_end
            inc_end = aend - cov_start + 1

    elif mode == 'span':
        if transcripts != None:
            tx = transcripts[anchor_id]

            tstart = map_transcript_start(tx, rstart, rend)
            if tstart == None:
                inc_start = None
                inc_end = None
            else:
                tend = tstart + rend - rstart

                if tx.strand == '-':
                    tstart_rev = transcript_lengths[anchor_id] - tend
                    tend_rev = transcript_lengths[anchor_id] - tstart
                    tstart = tstart_rev
                    tend = tend_rev

                tstart_pct = tstart/float(transcript_lengths[anchor_id])
                tend_pct = tend/float(transcript_lengths[anchor_id])

                inc_start = int(bins*tstart_pct)
                inc_end = int(0.5 + bins*tend_pct)
        else:
            cov_start = max(rstart, astart)
            cov_end = min(rend, aend)
            alength = aend - astart + 1

            if astrand == '+':                            
                cov_start_pct = (cov_start - astart) / float(alength)
                cov_end_pct = (cov_end - astart + 1) / float(alength)
            else:
                cov_start_pct = (aend - cov_end) / float(alength)
                cov_end_pct = (aend - cov_start + 1) / float(alength)

            inc_start = int(bins*cov_start_pct)
            inc_end = int(0.5 + bins*cov_end_pct)

    else:
        print >> sys.stderr, 'Unknown mode %s' % mode
        exit(1)

    return inc_start, inc_end


################################################################################
# initialize_coverage
################################################################################
def initialize_coverage(anchor_gff, mode, anchor_is_gtf, bins):
    print >> sys.stderr, 'Initializing coverage using anchor gff %s' % anchor_gff

    coverage = {}
    for line in open(anchor_gff):
        a = line.split('\t')

        chrom = a[0]
        start = int(a[3])
        end = int(a[4])
            
        if anchor_is_gtf:
            anchor_id = gff.gtf_kv(a[8])['transcript_id']
        else:
            anchor_id = '%s:%d-%d' % (chrom, start, end)
            
        if not anchor_id in coverage:
            if mode == 'span':
                coverage[anchor_id] = [0]*bins
            elif mode == 'mid':
                coverage[anchor_id] = [0]*(end-start+1)
            else:
                print >> sys.stderr, 'Unknown mode %s' % mode
                exit(1)

    print >> sys.stderr, '%d anchors found.' % len(coverage)

    return coverage


################################################################################
# map_transcript_start
#
# Given a transcript and read alignment, find the read's start (i.e. left) index
# on the transcript.
#
# Input
#  tx:     Transcript Gene object.
#  rstart: Read start (i.e. left) index on the genome.
#  rend:   Read end (i.e. right) index on the genome.
#
# Output
#  tstart: Read start (i.e. left) index on the transcript.
################################################################################
def map_transcript_start(tx, rstart, rend):
    tstart = 0

    for exon in tx.exons:
        # read is before exon
        if rend < exon.start:
            # done
            break

        # read is after exon
        elif exon.end < rstart:
            # add full exon length
            tstart += exon.end-exon.start+1

        # read is inside exon
        elif exon.start <= rstart and rend <= exon.end:
            # add partial exon length, and done
            tstart += rstart - exon.start
            break

        # read overlaps an exon boundary
        else:
            # read is incompatible with isoform so scrap it
            tstart = None
            break

    return tstart


################################################################################
# preprocess_anchors
#
# Input
#  anchor_gff:      Anchor GFF filename.
#  mode:            mid or span.
#  max_anchors:     Maximum # of anchors to use.
#  anchor_is_gtf:   True iff anchor_gff is a GTF file.
#  min_length:      Minimum anchor length.
#  window:          Window size around the midpoint.
#
# Output
#  prep_anchor_fd:  Preprocessed GFF file descriptor.
#  prep_anchor_gff: Preprocessed GFF filename.
################################################################################
def preprocess_anchors(anchor_gff, mode, max_anchors, anchor_is_gtf, min_length, window):
    # get lengths
    anchor_lengths = {}
    for line in open(anchor_gff):
        a = line.split('\t')

        if anchor_is_gtf:
            anchor_id = gff.gtf_kv(a[8])['transcript_id']
        else:
            anchor_id = (a[0], int(a[3]), int(a[4]))

        anchor_lengths[anchor_id] = anchor_lengths.get(anchor_id,0) + int(a[4])-int(a[3])+1

    # filter small
    if min_length != None:
        for anchor_id in anchor_lengths.keys():
            if anchor_lengths[anchor_id] < min_length:
                del anchor_lengths[anchor_id]

    # sample
    if max_anchors < len(anchor_lengths):
        anchors_chosen = set(random.sample(anchor_lengths.keys(), max_anchors))
    else:
        anchors_chosen = set(anchor_lengths.keys())

    # make new GFF
    prep_anchor_fd, prep_anchor_gff = tempfile.mkstemp()
    print >> sys.stderr, 'Opening tempfile %s for preprocessed anchors.' % prep_anchor_gff
    prep_anchor_out = open(prep_anchor_gff, 'w')

    for line in open(anchor_gff):
        a = line.split('\t')

        if anchor_is_gtf:
            anchor_id = gff.gtf_kv(a[8])['transcript_id']
        else:
            anchor_id = (a[0], int(a[3]), int(a[4]))

        if anchor_id in anchors_chosen:
            if mode == 'span':
                print >> prep_anchor_out, line,
            elif mode == 'mid':
                # standardize size
                start = int(a[3])
                end = int(a[4])
                mid = start + (end-start)/2
                a[3] = str(mid - window/2)
                a[4] = str(mid + window/2)
                a[-1] = a[-1].rstrip()

                if int(a[3]) > 0:
                    print >> prep_anchor_out, '\t'.join(a)
            else:
                print >> sys.stderr, 'Unknown mode %s' % mode
                exit(1)

    prep_anchor_out.close()

    return prep_anchor_fd, prep_anchor_gff


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
