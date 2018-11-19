#!/usr/bin/env python
from optparse import OptionParser
import subprocess

'''
make_ref_ml.py

Make machine learning friendly genome files, removing unplaced contigs,
and chrY.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fasta_file> <genome_file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide FASTA and genome files')
    else:
        fasta_file = args[0]
        genome_file = args[1]

    fasta_ml_file = fasta_file.replace('.fa', '.ml.fa')
    fasta_ml_out = open(fasta_ml_file, 'w')

    for line in open(fasta_file):
        if line[0] == '>':
            keep_chr = True
            header = line[1:]
            keep_chr = filter_chr(header)
        if keep_chr:
            print(line, file=fasta_ml_out, end='')

    fasta_ml_out.close()

    subprocess.call('samtools faidx %s' % fasta_ml_file, shell=True)


    genome_ml_file = genome_file.replace('.genome', '.ml.genome')
    genome_ml_file = open(genome_ml_file, 'w')

    for line in open(genome_file):
        header = line.split()[0]
        keep_chr = filter_chr(header)
        if keep_chr:
            print(line, file=genome_ml_file, end='')

    genome_ml_file.close()


def filter_chr(header):
    keep_chr = True
    if header.find('chrUn') != -1:
        keep_chr = False
    elif header.find('random') != -1:
        keep_chr = False
    elif header.find('hap') != -1:
        keep_chr = False
    elif header.find('alt') != -1:
        keep_chr = False
    elif header.find('KI270') != -1:
        keep_chr = False
    elif header.find('GL000') != -1:
        keep_chr = False
    elif header.find('JH584') != -1:
        keep_chr = False
    elif header.find('GL456') != -1:
        keep_chr = False
    elif header.rstrip() in ['chrM','chrMT','chrY']:
        keep_chr = False
    return keep_chr

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
