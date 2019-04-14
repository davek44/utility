#!/usr/bin/env python
from optparse import OptionParser
import pdb
import os

import pandas as pd

from basenji.emerald import EmeraldVCF

'''
vcf_ld.py

Transform an input VCF to add all linked variants above some threshold.
Makes use of Emerald for LD queries.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_vcf> <out_vcf>'
    parser = OptionParser(usage)
    parser.add_option('-l','--ld', dest='ld_t',
        default=0.8, type='float',
        help='LD threshold to include SNP [Default: %default]')
    parser.add_option('-r', dest='refpanel_stem',
        default='%s/popgen/1000G/phase3/eur/1000G.EUR.QC' % os.environ['HG19'],
        help='Reference panel chromosome VCF stem [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input and output VCF files')
    else:
        in_vcf_file = args[0]
        out_vcf_file = args[1]

    # initialize reference panel
    refp_em = EmeraldVCF(options.refpanel_stem)

    # retrieve all SNPs in LD
    all_snps_df = []

    # initialize VCFs
    in_vcf_open = open(in_vcf_file)
    out_vcf_open = open(out_vcf_file, 'w')

    # hash SNPs by chromosome
    for line in in_vcf_open:
        if line[0] == '#':
            # print header
            print(line, end='', file=out_vcf_open)

        else:
            a = line.split()
            chrm = a[0]
            pos = int(a[1])
            rsid = a[2]

            # query LD SNPs
            snp_df = refp_em.query_ld(rsid, chrm, pos,
                                      options.ld_t, return_pos=True)

            if snp_df.shape[0] == 0:
                print('WARNING: %s not found in reference panel.' % rsid)
            else:
                # set SNP id index
                snp_df.set_index('snp', inplace=True)

                # fetch VCF lines
                pos_start = snp_df.pos.iloc[0]
                pos_end = snp_df.pos.iloc[-1]
                for snp_rec in refp_em.fetch(chrm, pos_start-1, pos_end):
                    if snp_rec.id in snp_df.index:
                        snp_str = snp_rec.__str__().rstrip()
                        snp_str += '=%s;LD=%.2f' % (rsid, snp_df.loc[snp_rec.id].r)
                        print(snp_str, file=out_vcf_open)

    out_vcf_open.close()
    in_vcf_open.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
