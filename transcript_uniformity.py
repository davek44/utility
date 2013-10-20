#!/usr/bin/env python
from optparse import OptionParser
import math, os, pdb, subprocess, sys, tempfile
import pysam
import gff, stats
import clip_peaks

################################################################################
# transcript_uniformity.py
#
# Measure the uniformity of read coverage using the index of dispersion on the
# transcripts defined by a given GTF file.
################################################################################

# globals that are a pain to pass around
clip_peaks.out_dir = None

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam> <ref_gtf>'
    parser = OptionParser(usage)

    # IO options
    parser.add_option('-o', dest='out_dir', default='uniform', help='Output directory [Default: %default]')

    # window options
    parser.add_option('-w', dest='window_size', type='int', default=25, help='Window size for counting [Default: %default]')
    parser.add_option('-i', '--ignore', dest='ignore_gff', help='Ignore reads overlapping overlapping troublesome regions in the given GFF file')
    parser.add_option('-u', '--unstranded', dest='unstranded', action='store_true', default=False, help='Sequencing is unstranded [Default: %default]')

    # cufflinks options
    parser.add_option('--cuff_done', dest='cuff_done', action='store_true', default=False, help='The Cufflinks run to estimate the model parameters is already done [Default: %default]')
    parser.add_option('-t', dest='threads', type='int', default=2, help='Number of threads to use [Default: %default]')

    # debug options
    parser.add_option('-v', '--verbose', dest='verbose', action='store_true', default=False, help='Verbose output [Default: %default]')
    parser.add_option('-g', '--gene', dest='gene_only', help='Call peaks on the specified gene only')
    #parser.add_option('--print_windows', dest='print_windows', default=False, action='store_true', help='Print statistics for all windows [Default: %default]')

    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error(usage)
    else:
        bam = args[0]
        ref_gtf = args[1]

    clip_peaks.out_dir = options.out_dir

    if not os.path.isdir(clip_peaks.out_dir):
        os.mkdir(clip_peaks.out_dir)

    ############################################
    # parameterize
    ############################################
    if not options.cuff_done:
        # make a new gtf w/ unspliced RNAs
        update_ref_gtf = clip_peaks.prerna_gtf(ref_gtf)

        subprocess.call('cufflinks -o %s -p %d -G %s %s' % (clip_peaks.out_dir, options.threads, update_ref_gtf, bam), shell=True)

    # store transcripts
    transcripts = clip_peaks.read_genes('%s/transcripts.gtf'%clip_peaks.out_dir, key_id='transcript_id')

    # merge overlapping genes
    g2t_merge, antisense_clusters = clip_peaks.merged_g2t('%s/transcripts.gtf'%clip_peaks.out_dir, options.unstranded)

    if options.unstranded:
        # alter strands
        clip_peaks.ambiguate_strands(transcripts, g2t_merge, antisense_clusters)

    # set transcript FPKMs
    clip_peaks.set_transcript_fpkms(transcripts, clip_peaks.out_dir, missing_fpkm=0)

    # possibly limit genes to examine
    if options.gene_only:
        gene_ids = []
        for gids in g2t_merge.keys():
            if options.gene_only in gids.split(','):
                gene_ids.append(gids)
        if len(gene_ids) == 0:
            print >> sys.stderr, 'gene_id %s not found' % options.gene_only
            exit(1)
    else:
        gene_ids = g2t_merge.keys()


    ############################################
    # filter BAM
    ############################################
    if options.ignore_gff:
        bam_ignore_fd, bam_ignore_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
        subprocess.call('intersectBed -v -abam %s -b %s > %s' % (bam, options.ignore_gff, bam_ignore_file), shell=True)
        bam = bam_ignore_file

    ############################################
    # process genes
    ############################################
    # index
    subprocess.call('samtools index %s' % bam, shell=True)

    # initialize stats
    table_out = open('%s/uniformity_table.txt' % clip_peaks.out_dir, 'w')
    id_list = []
    fpkm_list = []

    # open bam
    bam_in = pysam.Samfile(bam, 'rb')
    
    # for each gene
    for gene_id in gene_ids:
        # make a more focused transcript hash for this gene
        gene_transcripts = {}
        for tid in g2t_merge[gene_id]:
            gene_transcripts[tid] = transcripts[tid]

        # obtain basic gene attributes
        (gchrom, gstrand, gstart, gend) = clip_peaks.gene_attrs(gene_transcripts)

        # initialize window counts
        transcript_isoform_counts = {}
        for tid in gene_transcripts:
            transcript_isoform_counts[tid] = []

        # choose a single event position and weight the reads
        read_pos_weights = clip_peaks.position_reads(bam_in, gchrom, gstart, gend, gstrand, mapq_zero=True)

        # process read alignments
        for (pos, weight, mm) in read_pos_weights:
            # map pos to isoforms
            iso_pos = {}
            for tid in gene_transcripts:
                iso_pos[tid] = isoform_position(gene_transcripts[tid], pos)

            # sum fpkms for hit isoforms
            fpkm_sum = sum([gene_transcripts[tid].fpkm for tid in gene_transcripts if iso_pos[tid] != None])

            if fpkm_sum <= 0:
                pass
                #print >> sys.stderr, 'No FPKM for %s at %d' % (gene_id,pos)
            else:
                # distribute read to isoform counts
                for tid in gene_transcripts:
                    if iso_pos[tid] != None:
                        win_i = int(iso_pos[tid] / options.window_size)
                        while win_i >= len(transcript_isoform_counts[tid]):
                            transcript_isoform_counts[tid].append(0)
                        transcript_isoform_counts[tid][win_i] += weight*gene_transcripts[tid].fpkm/fpkm_sum

        # compute window stats
        for tid in gene_transcripts:
            if gene_transcripts[tid].fpkm > 1 and len(transcript_isoform_counts[tid]) > 5:
                u, sd = stats.mean_sd(transcript_isoform_counts[tid][:-1])
                if u > 0:
                    id_list.append(sd*sd/u)
                    fpkm_list.append(gene_transcripts[tid].fpkm)

                    cols = (tid, gene_transcripts[tid].fpkm, len(transcript_isoform_counts[tid])-1, u, sd, id_list[-1])
                    print >> table_out, '%-20s  %8.2f  %6d  %7.2f  %7.2f  %5.3f' % cols        

    bam_in.close()
    table_out.close()

    ############################################
    # summary stats
    ############################################
    median = stats.median(id_list)
    mean = stats.mean(id_list)

    fpkm_cv_sum = sum([id_list[i]*fpkm_list[i] for i in range(len(id_list))])
    fpkm_sum = sum(fpkm_list)
    fpkm_mean = fpkm_cv_sum / fpkm_sum

    logfpkm_cv_sum = sum([id_list[i]*math.log(fpkm_list[i]+1,2) for i in range(len(id_list))])
    logfpkm_sum = sum([math.log(f+1,2) for f in fpkm_list])
    logfpkm_mean = logfpkm_cv_sum / logfpkm_sum

    # print
    print 'Median:                %7.4f' % median
    print 'Mean:                  %7.4f' % mean
    print 'FPKM-weighted mean:    %7.4f' % fpkm_mean
    print 'logFPKM-weighted mean: %7.4f' % logfpkm_mean

    # clean cufflinks output
    if not options.cuff_done:
        os.remove(update_ref_gtf)
        os.remove('%s/skipped.gtf' % clip_peaks.out_dir)
        os.remove('%s/genes.fpkm_tracking' % clip_peaks.out_dir)

    if options.ignore_gff:
        os.close(bam_ignore_fd)
        os.remove(bam_ignore_file)


################################################################################
# isoform_position
#
# Map the given genomic position to a relative isoform position.
#
# Input
#  transcript: Gene class object representing an individual isoform.
#  genome_pos: Genomic position.
#
# Output:
#  iso_pos:    Isoform position.
################################################################################
def isoform_position(transcript, genome_pos):
    iso_pos = 0
    for exon_i in range(len(transcript.exons)):
        if genome_pos < transcript.exons[exon_i].start:
            return None
        elif genome_pos < transcript.exons[exon_i].end:
            return iso_pos + genome_pos - transcript.exons[exon_i].start
        else:
            iso_pos += transcript.exons[exon_i].end - transcript.exons[exon_i].start + 1
    return None


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
