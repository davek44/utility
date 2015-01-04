#!/usr/bin/env python
from optparse import OptionParser
from numpy import mean
import os, pdb, shutil, subprocess
import cuffdiff, fdr, gff, ggplot, math, stats, te

################################################################################
# te_diff.py
#
# Compute stats and plot differential expression fold changes for genes
# w/ and w/o each TE family.
#
# Written to operate on transcript_id's, assuming ref_gtf has been filtered
# to have a single isoform per gene.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf> <diff>'
    parser = OptionParser(usage)
    parser.add_option('-m', dest='max_stat', default=None, type='float', help='Maximum stat for plotting [Default: %default]')
    parser.add_option('-o', dest='out_dir', default='te_diff', help='Output directory [Default: %default]')
    parser.add_option('-c', dest='scale', default=1, type='float', help='CDF plot scale [Default: %default]')
    parser.add_option('-t', dest='te_gff', default='%s/hg19.fa.out.tp.gff'%os.environ['MASK'])

    parser.add_option('-s', dest='spread_factor', default=None, type='float', help='Allow multiplicative factor between the shortest and longest transcripts, used to filter [Default: %default]')
    parser.add_option('-l', dest='spread_lower', default=None, type='float', help='Allow multiplicative factor between median and shortest transcripts [Defafult: %default]')
    parser.add_option('-u', dest='spread_upper', default=None, type='float', help='Allow multiplicative factor between median and longest transcripts [Defafult: %default]')

    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide .gtf and .diff files')
    else:
        ref_gtf = args[0]
        diff_file = args[1]

    # clean plot directory
    if os.path.isdir(options.out_dir):
        shutil.rmtree(options.out_dir)
    os.mkdir(options.out_dir)

    if options.spread_factor or options.spread_lower or options.spread_upper:
        # filter for similar length

        if options.spread_factor:
            options.spread_lower = math.sqrt(options.spread_factor)
            options.spread_upper = math.sqrt(options.spread_factor)

        spread_gtf = '%s/spread_factor.gtf' % options.out_dir
        gff.length_filter(ref_gtf, spread_gtf, options.spread_lower, options.spread_lower)

        ref_gtf = spread_gtf

    ##################################################
    # hash TEs -> genes
    ##################################################
    te_genes = te.hash_repeats_genes(ref_gtf, options.te_gff, gene_key='transcript_id', add_star=True, stranded=True)

    ##################################################
    # hash genes -> RIP diff
    ##################################################
    gene_diff = cuffdiff.hash_diff(diff_file, stat='fold', max_stat=options.max_stat, sample_first='input')

    ##################################################
    # compute stats and make plots
    ##################################################
    table_lines, pvals = compute_stats(te_genes, gene_diff, ref_gtf, options.out_dir, options.scale)

    # perform multiple hypothesis correction
    qvals = fdr.ben_hoch(pvals)

    table_out = open('%s/table.txt' % options.out_dir, 'w')
    for i in range(len(table_lines)):
        print >> table_out, '%s %10.2e' % (table_lines[i],qvals[i])
    table_out.close()

    if options.spread_factor or options.spread_lower or options.spread_upper:
        os.remove(spread_gtf)


################################################################################
# cdf_plot
#
# Input
#  te_key:
#  te_diffs:
#  note_diffs:
#  out_pdf:
################################################################################
def cdf_plot(te_key, te_diffs, note_diffs, out_pdf, scale):
    rep, fam, orient = te_key

    # name plot
    if fam == '-':
        label = 'dTE-RNAs/%s' % orient
    elif fam == '*':
        label = 'TE-RNAs/%s' % orient
    elif rep == '*':
        label = '%s-RNAs/%s' % (fam,orient)
    else:
        label = '%s-RNAs/%s' % (rep,orient)

    # construct data frame
    df = {}
    df['diff'] = note_diffs + te_diffs
    df['class'] = ['d%s' % label]*len(note_diffs) + [label]*len(te_diffs)

    ggplot.plot('%s/te_diff.r' % os.environ['RDIR'], df, [out_pdf, scale])


################################################################################
# compute_stats
#
# Input
#  te_genes:
#  gene_diff:
#  ref_gtf:
#  plot_dir:
#
# Output
#  table_lines:
#  pvals:
################################################################################
def compute_stats(te_genes, gene_diff, ref_gtf, plot_dir, scale):
    # focus on GTF genes
    gtf_genes = set()
    for line in open(ref_gtf):
        a = line.split('\t')
        tid = gff.gtf_kv(a[8])['transcript_id']
        gtf_genes.add(tid)    

    pvals = []
    table_lines = []

    for sample_key in gene_diff:
        sample1, sample2 = sample_key

        stat_genes = list(gtf_genes & set(gene_diff[sample_key]))

        for te_key in te_genes:
            repeat, family, orient = te_key

            te_diffs = [gene_diff[sample_key][tid] for tid in stat_genes if tid in te_genes[te_key]]
            if len(te_diffs) > 0:
                note_diffs = [gene_diff[sample_key][tid] for tid in stat_genes if tid not in te_genes[te_key]]

                te_mean = mean(te_diffs)
                note_mean = mean(note_diffs)

                if len(te_diffs) > 5:
                    z, p = stats.mannwhitneyu(te_diffs, note_diffs)
                else:
                    z = 0
                    p = 1

                pvals.append(p)

                cols = (repeat, family, orient, sample1, sample2, len(te_diffs), te_mean, len(note_diffs), note_mean, z, p)
                table_lines.append('%-17s  %-17s  %1s  %-10s  %-10s  %6d  %9.2f  %6d  %9.2f  %8.2f  %10.2e' % cols)

                # plot ...
                if repeat in ['*'] and family in ['*','LINE/L1','SINE/Alu','LTR/ERV1','LTR/ERVL-MaLR','LINE/L2','LTR/ERVL','SINE/MIR','DNA/hAT-Charlie','LTR/ERVK','DNA/TcMar-Tigger']:
                    repeat_plot = repeat.replace('/','-').replace('*','X')
                    family_plot = family.replace('/','-').replace('*','X')
                    out_pdf = '%s/%s_%s_%s_%s-%s.pdf' % (plot_dir, repeat_plot, family_plot, orient, sample1, sample2)
                    cdf_plot(te_key, te_diffs, note_diffs, out_pdf, scale)

    return table_lines, pvals


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
