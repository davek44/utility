#!/usr/bin/env python
from optparse import OptionParser

'''
bim_vcf.py

Convert variants in a Plink .bim file to .vcf
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_bim> <out_vcf>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input BIM and output VCF.')
    else:
        in_bim_file = args[0]
        out_vcf_file = args[1]

    # open out VCF
    out_vcf_open = open(out_vcf_file, 'w')

    # print header
    print('##fileformat=VCFv4.2', file=out_vcf_open)
    cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    print('\t'.join(cols), file=out_vcf_open)

    # parse BIM
    for line in open(in_bim_file):
        a = line.split()
        chrom = a[0]
        snp_id = a[1]
        pos = a[3]
        a1 = a[4]
        a2 = a[5]

        cols = [chrom, pos, snp_id, a2, a1, '.', '.', '.']
        print('\t'.join(cols), file=out_vcf_open)

    out_vcf_open.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
