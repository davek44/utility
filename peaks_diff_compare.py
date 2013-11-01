#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess
import ggplot, gff, stats

################################################################################
# peaks_diff_compare.py
#
# Compare RNAs bound according to some peak calls to a cuffdiff run.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <peaks gff> <diff>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='control_fpkm_file', help='Control FPKM tracking file')
    parser.add_option('-g', dest='ref_gtf', default='%s/gencode.v18.annotation.gtf'%os.environ['GENCODE'])
    parser.add_option('-o', dest='output_pre', default='', help='Output prefix [Default: %default]')
    parser.add_option('--sample1', dest='sample1', help='Sample_1 name in cuffdiff')
    parser.add_option('--sample2', dest='sample2', help='Sample_2 name in cuffdiff')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide peaks GFF and .diff file')
    else:
        peaks_gff = args[0]
        diff_file = args[1]

    # find expressed genes in peak calls
    silent_genes = set()
    if options.control_fpkm_file:
        silent_genes = find_silent(options.control_fpkm_file)

    # find peak bound genes
    peak_genes = set()
    p = subprocess.Popen('intersectBed -s -u -a %s -b %s' % (options.ref_gtf,peaks_gff), shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        peak_genes.add(gff.gtf_kv(line.split('\t')[8])['gene_id'])
    p.communicate()

    # store test stats
    bound_tstats = []
    unbound_tstats = []

    diff_in = open(diff_file)
    line = diff_in.readline()
    for line in diff_in:
        a = line.split('\t')

        gene_id = a[0]
        sample1 = a[4]
        sample2 = a[5]
        status = a[6]
        fpkm1 = float(a[7])
        fpkm2 = float(a[8])
        tstat = float(a[10])

        if status == 'OK' and not math.isnan(tstat):
            if options.sample1 in [None,sample1] and options.sample2 in [None,sample2]:
                if gene_id in peak_genes:
                    bound_tstats.append(tstat)
                else:
                    if not gene_id in silent_genes:
                        unbound_tstats.append(tstat)

    print '%d silent genes' % len(silent_genes)
    print '%d bound genes' % len(bound_tstats)
    print '%d unbound genes' % len(unbound_tstats)

    # perform statistical test
    z, p = stats.mannwhitneyu(bound_tstats, unbound_tstats)
    print z, p

    # construct data frame
    df_dict = {'Peak':(['Yes']*len(bound_tstats) + ['No']*len(unbound_tstats)),
               'Test_stat':bound_tstats+unbound_tstats}

    r_script = '%s/peaks_diff_compare.r' % os.environ['GGPLOT']

    output_pdf = '%s_dens.pdf' % options.output_pre

    ggplot.plot(r_script, df_dict, [output_pdf])
    

################################################################################
# find_silent
#
# Input:
#  control_fpkm_file: Cufflinks FPKM file.
#  silent_fpkm:       FPKM threshold to call a gene silent.
#
# Output:
#  silent_genes:      Set of silent gene_id's.
################################################################################
def find_silent(control_fpkm_file, silent_fpkm=0.5):
    # get fpkms (possibly from an isoform file)
    gene_fpkms = {}
    control_fpkm_in = open(control_fpkm_file)
    control_fpkm_in.readline()
    for line in control_fpkm_in:
        a = line.split('\t')
        gene_id = a[3]
        fpkm = float(a[9])
        gene_fpkms[gene_id] = gene_fpkms.get(gene_id,0) + fpkm
    control_fpkm_in.close()

    silent_genes = set()
    for gene_id in gene_fpkms:
        if gene_fpkms[gene_id] < silent_fpkm:
            silent_genes.add(gene_id)

    return silent_genes

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
