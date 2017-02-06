#!/usr/bin/env python
from optparse import OptionParser
import gc
import os

import numpy as np
import pandas as pd

################################################################################
# vcf_refpanel.py
#
# Filter a sorted VCF file for only SNPs present and matching the 1kG reference
# panel. Flip SNPs with ref and alt reversed.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_vcf_file> <out_vcf_file>'
    parser = OptionParser(usage)
    parser.add_option('-r', dest='refpanel_prefix', default='%s/popgen/1000G_EUR_Phase3_plink/1000G.EUR.QC' % os.environ['HG19'], help='Reference panel directory with Plink files by chromosome.')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input and output VCF files')
    else:
        in_vcf_file = args[0]
        out_vcf_file = args[1]

    pd.options.mode.chained_assignment = None

    # read VCF
    try:
        df_vcf = pd.read_table(in_vcf_file, skiprows=1, names=['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter'])
    except:
        df_vcf = pd.read_table(in_vcf_file, skiprows=1, names=['chr', 'pos', 'id', 'ref', 'alt'])

    # open output VCF
    vcf_out = open(out_vcf_file, 'w')

    # print header
    vcf_in = open(in_vcf_file)
    line = vcf_in.readline()
    while line.startswith('#'):
        print(line, file=vcf_out, end='')
        line = vcf_in.readline()
    vcf_in.close()

    chroms = set(df_vcf.chr)
    for chrom in chroms:
        # clean up past views
        gc.collect()

        # read reference chromosome table
        chrom_num = chrom[3:]
        bim_chrom = '%s.%s.bim' % (options.refpanel_prefix, chrom_num)
        if not os.path.isfile(bim_chrom):
            print('Skipping %s' % chrom)
            continue
        df_ref_chrom = pd.read_table(bim_chrom, names=['chr', 'id', 'z', 'pos', 'ref', 'alt'])

        # filter VCF table for chromosome
        df_vcf_chrom = df_vcf.loc[df_vcf.chr == chrom]

        # determine shared positions
        positions_vcf_shared = np.in1d(df_vcf_chrom.pos, df_ref_chrom.pos)
        df_vcf_chrom_shared = df_vcf_chrom.loc[positions_vcf_shared]
        positions_ref_shared = np.in1d(df_ref_chrom.pos, df_vcf_chrom.pos)
        df_ref_chrom_shared = df_ref_chrom.loc[positions_ref_shared]
        
        # determine matching alleles
        match_ref = (np.array(df_vcf_chrom_shared.ref) == np.array(df_ref_chrom_shared.ref))
        match_alt = (np.array(df_vcf_chrom_shared.alt) == np.array(df_ref_chrom_shared.alt))
        match_asis = match_ref & match_alt

        # determine flip alleles
        match_ref_alt = (np.array(df_vcf_chrom_shared.ref) == np.array(df_ref_chrom_shared.alt))
        match_alt_ref = (np.array(df_vcf_chrom_shared.alt) == np.array(df_ref_chrom_shared.ref))
        match_flip = match_ref_alt & match_alt_ref

        # flip
        tmp_allele = df_vcf_chrom_shared.ref[match_flip]
        df_vcf_chrom_shared.ref[match_flip] = df_vcf_chrom_shared.alt[match_flip]
        df_vcf_chrom_shared.alt[match_flip] = tmp_allele

        # print
        match_valid = match_asis | match_flip
        df_vcf_chrom_shared.loc[match_valid].to_csv(vcf_out, sep='\t', index=False, header=False)

    # close output VCF
    vcf_out.close()


def print_snp(snp):
    cols = [snp.chr, str(snp.pos), snp.id, snp.ref, snp.alt]
    if 'qual' in snp:
        cols += [snp.qual, snp['filter']]
    print('\t'.join(cols))


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
