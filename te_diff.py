#!/usr/bin/env python
from optparse import OptionParser
from collections import Counter
import os, pdb, shutil, subprocess
import fdr, gff, ggplot, math, stats, te

################################################################################
# te_cuffdiff.py
#
# Compute stats and plot differential expression fold changes for genes
# w/ and w/o each TE family.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf> <diff>'
    parser = OptionParser(usage)
    parser.add_option('-o', dest='out_dir', default='te_diff', help='Output directory [Default: %default]')
    parser.add_option('-t', dest='te_gff', default='%s/hg19.fa.out.tpf.gff'%os.environ['MASK'])
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide .gtf and .diff files')
    else:
        gtf_file = args[0]
        diff_file = args[1]

    # hash genes -> TEs
    gene_tes = te.hash_genes_repeats(gtf_file, options.te_gff, gene_key='transcript_id', add_star=True, stranded=True)

    # create a fake family for unrepetitive genes
    for line in open(gtf_file):
        a = line.split('\t')
        gene_id = gff.gtf_kv(a[8])['transcript_id']
        if not gene_id in gene_tes:
            gene_tes[gene_id] = set([('-','-','*')])

    # get diffs stats
    gene_diffs, te_diffs = get_diff_stats(diff_file, gene_tes)

    # clean plot directory
    if os.path.isdir(options.out_dir):
        shutil.rmtree(options.out_dir)
    os.mkdir(options.out_dir)

    # stats
    table_lines, pvals = compute_stats(te_diffs, gene_diffs, options.out_dir)

    # perform multiple hypothesis correction
    qvals = fdr.ben_hoch(pvals)

    table_out = open('%s/table.txt' % options.out_dir, 'w')
    for i in range(len(table_lines)):
        print >> table_out, '%s %10.2e' % (table_lines[i],qvals[i])
    table_out.close()


################################################################################
# cdf_plot
################################################################################
def cdf_plot(te_or, w_te, wo_te, out_pdf):
    rep, fam, orient = te_or

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
    df['fold'] = wo_te + w_te
    df['class'] = ['d%s' % label]*len(wo_te) + [label]*len(w_te)

    ggplot.plot('te_diff.r', df, [out_pdf])


################################################################################
# compute_stats
################################################################################
def compute_stats(te_diffs, gene_diffs, plot_dir):
    pvals = []
    table_lines = []

    for te_or in te_diffs:
        rep, fam, orient = te_or
        
        for sample_key in te_diffs[te_or]:        
            sample1, sample2 = sample_key

            # if enough data
            if len(te_diffs[te_or][sample_key]) >= 10:
                wo_te = list((gene_diffs[sample_key] - te_diffs[te_or][sample_key]).elements())
                w_te = list(te_diffs[te_or][sample_key].elements())

                wo_mean = stats.mean(wo_te)
                w_mean = stats.mean(w_te)

                z, p = stats.mannwhitneyu(w_te, wo_te)

                cols = (rep, fam, orient, sample1, sample2, len(w_te), w_mean, wo_mean, z, p)
                table_lines.append('%-17s %-17s  %1s  %-10s %-10s %6d %9.2f %9.2f %8.2f %10.2e' % cols)

                pvals.append(p)

                # plot ...
                if rep in ['*'] and fam in ['*','LINE/L1','SINE/Alu','LTR/ERV1','LTR/ERVL-MaLR','LINE/L2','LTR/ERVL','SINE/MIR','DNA/hAT-Charlie','LTR/ERVK','DNA/TcMar-Tigger']:
                    out_pdf = '%s/%s_%s_%s_%s-%s.pdf' % (plot_dir,rep.replace('/','-'),fam.replace('/','-'),orient,sample1,sample2)
                    cdf_plot(te_or, w_te, wo_te, out_pdf)

    return table_lines, pvals


################################################################################
# get_diff_stats
################################################################################
def get_diff_stats(diff_file, gene_tes):
    # initialize diff counters
    gene_diffs = {}
    te_diffs = {}

    # read diff file
    diff_in = open(diff_file)
    headers = diff_in.readline()
    line = diff_in.readline()
    while line:
        a = line.split('\t')

        gene_id = a[0]
        sample1 = a[4]
        sample2 = a[5]
        status = a[6]
        fpkm1 = float(a[7])
        fpkm2 = float(a[8])
        fold = float(a[9])
        tstat = float(a[10])
        sig = a[-1].rstrip()

        if sample2 == 'input':
            sample1, sample2 = sample2, sample1
            fpkm1, fpkm2 = fpkm2, fpkm1
            fold *= -1
            tstat *= -1

        # cap fold/tstat
        fold = min(fold, 6)
        fold = max(fold, -6)
        tstat = min(tstat, 6)
        tstat = max(tstat, -6)

        if gene_id in gene_tes and status == 'OK' and not math.isnan(tstat):
            # save for global
            #gene_diffs.setdefault((sample1,sample2),Counter())[tstat] += 1
            gene_diffs.setdefault((sample1,sample2),Counter())[fold] += 1

            # save for TEs
            for te_or in gene_tes[gene_id]:
                if not te_or in te_diffs:
                    te_diffs[te_or] = {}
                #te_diffs[te_or].setdefault((sample1,sample2),Counter())[tstat] += 1
                te_diffs[te_or].setdefault((sample1,sample2),Counter())[fold] += 1

        line = diff_in.readline()
    diff_in.close()

    return gene_diffs, te_diffs


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
