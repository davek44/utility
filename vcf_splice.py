#!/usr/bin/env python
from optparse import OptionParser
import os
import pdb

import pybedtools

'''
vcf_splice.py

Add splice site distance INFO column to a VCF file.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <in_vcf_file> <out_vcf_file>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='splice_gff_file', default='%s/genes/gencode28/gencode.v28.basic.annotation.splice.gff' % os.environ['HG38'])
    # parser.add_option('-g', dest='splice_gff_file', default='%s/genes/gencode28/gencode_basic_splice.gff' % os.environ['HG19'])
    parser.add_option('-t', dest='filter_t',
        default=None, type='int',
        help='Filter out variants less than the given distance threshold [Default: %default]')
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
            print('##FORMAT=<ID=SS,Number=1,Type=Integer,Description="Splice site distance">', file=out_vcf_open)
        print(line, end='', file=out_vcf_open)
        line = in_vcf_open.readline()
    in_vcf_open.close()

    # intersect
    in_vcf_bedtool = pybedtools.BedTool(in_vcf_file)
    splice_bedtool = pybedtools.BedTool(options.splice_gff_file)

    for closest_a in in_vcf_bedtool.closest(splice_bedtool, d=True, t='first'):
        a = closest_a[:8]
        splice_distance = int(closest_a[-1])
        if a[-1] == '.':
            a[-1] = 'SS=%s' % str(splice_distance)
        else:

            a[-1] += ';SS=%s' % str(splice_distance)

        if options.filter_t is None or splice_distance >= options.filter_t:
            print('\t'.join(a), file=out_vcf_open)

    # close
    out_vcf_open.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
