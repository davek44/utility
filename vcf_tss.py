#!/usr/bin/env python
from optparse import OptionParser
import os
import pdb

import pybedtools

'''
vcf_tss.py

Add TSS distance INFO column to a VCF file.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_vcf_file> <out_vcf_file>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='tss_gff_file', default='%s/genes/gencode28/gencode.v28.basic.annotation.tss.gff' % os.environ['HG38'])
    # parser.add_option('-g', dest='tss_gff_file', default='%s/genes/gencode28/gencode_basic_tss.gff' % os.environ['HG19'])
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input and output VCF files')
    else:
        in_vcf_file = args[0]
        out_vcf_file = args[1]

    # open files
    in_vcf_open = open(in_vcf_file)
    out_vcf_open = open(out_vcf_file, 'w')

    # print header
    line = in_vcf_open.readline()
    while line.startswith('#'):
        if line.startswith('#CHROM'):
            # add new INFO description first
            print('##FORMAT=<ID=TS,Number=1,Type=Integer,Description="TSS distance">', file=out_vcf_open)
        print(line, end='', file=out_vcf_open)
        line = in_vcf_open.readline()
    in_vcf_open.close()

    # intersect
    in_vcf_bedtool = pybedtools.BedTool(in_vcf_file)
    tss_bedtool = pybedtools.BedTool(options.tss_gff_file)

    for closest_a in in_vcf_bedtool.closest(tss_bedtool, d=True, t='first'):
        a = closest_a[:8]
        if a[-1] == '.':
            a[-1] = 'TS=%s' % closest_a[-1]
        else:
            a[-1] += ';TS=%s' % closest_a[-1]
        print('\t'.join(a), file=out_vcf_open)

    # close
    out_vcf_open.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
