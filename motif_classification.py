#!/usr/bin/env python
from optparse import OptionParser
import numpy as np
import os
import shutil
import subprocess

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

################################################################################
# motif_classification.py
#
# Classify sequences in a BED/FASTA file using a given motif database.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bed/fasta>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='genome', default='HG19', help='Genome [Default: %default]')
    parser.add_option('-l', dest='model', default='logistic', help='sklearn model [Default: %default]')
    parser.add_option('-m', dest='motifs_file', default='%s/motif_databases/CIS-BP/Homo_sapiens.meme'%os.environ['MEME'], help='Motif database [Default: %default]')
    parser.add_option('-n', dest='neg_file')
    parser.add_option('-o', dest='out_dir', default='motif_class', help='Output directory')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide BED/FASTA file')
    else:
        input_file = args[0]

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    #######################################################
    # prepare datasets
    #######################################################

    # prepare positive set
    input_fasta = prep_input(input_file, options.genome, options.out_dir)

    # prepare negative set
    if options.neg_file:
        neg_fasta = prep_input(options.neg_file, options.genome, options.out_dir)
    else:
        neg_fasta = shuffle_input(input_file, options.genome, options.out_dir)

    #######################################################
    # compute motif features
    #######################################################
    X_pos = map_motifs(input_fasta, options.motifs_file, '%s/fimo_pos'%options.out_dir)
    X_neg = map_motifs(neg_fasta, options.motifs_file, '%s/fimo_neg'%options.out_dir)

    X = np.vstack([X_pos, X_neg])
    Y = np.array([1]*X_pos.shape[0] + [0]*X_neg.shape[0])

    #######################################################
    # fit model parameters
    #######################################################
    model = LogisticRegression()
    model.fit(X,Y)

    # test
    preds = model.predict_proba(X)
    auc = roc_auc_score(Y, preds[:,1])
    print('AUC: %f' % auc)


    #######################################################
    # analyze model
    #######################################################
    print(model.coef_.shape)

    coefs_out = open('%s/motif_coefs.txt' % options.out_dir, 'w')

    mi = 0
    for line in open(options.motifs_file):
        if line.startswith('MOTIF'):
            a = line.split()
            motif_id = a[1]
            protein = trim_parens(a[2].split('_')[0])

            print('%s\t%s\t%f' % (motif_id, protein, model.coef_[0,mi]), file=coefs_out)

            mi += 1

    coefs_out.close()

    #######################################################
    # clean up
    #######################################################
    input_suffix = os.path.splitext(input_file)[1]
    if input_suffix not in ['fa','fasta']:
        os.remove(input_fasta)

    if options.neg_file is None or os.path.splitext(options.neg_file)[1] not in ['fa', 'fasta']:
        os.remove(neg_fasta)


def map_motifs(fasta_file, motifs_file, out_dir):
    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)

    # run fimo
    cmd = 'fimo --bgfile motif-file -o %s %s %s' % (out_dir, motifs_file, fasta_file)
    subprocess.call(cmd, shell=True)

    # hash FASTA to indexes
    sequence_indexes = {}
    si = 0
    for line in open(fasta_file):
        if line[0] == '>':
            seq_id = line[1:].rstrip()
            sequence_indexes[seq_id] = si
            si += 1

    # hash motifs to indexes
    motif_indexes = {}
    mi = 0
    for line in open(motifs_file):
        if line.startswith('MOTIF'):
            motif_id = line.split()[1]
            motif_indexes[motif_id] = mi
            mi += 1

    # construct array from fimo output
    X = np.zeros((len(sequence_indexes), len(motif_indexes)))

    for line in open('%s/fimo.txt' % out_dir):
        if not line.startswith('#'):
            a = line.split()
            motif_id, seq_id = line.split()[:2]
            si = sequence_indexes[seq_id]
            mi = motif_indexes[motif_id]
            X[si,mi] += 1

    return X


def prep_input(input_file, genome, out_dir):
    ''' prepare input file as fasta '''

    input_suffix = os.path.splitext(input_file)[1]
    if input_suffix in ['fa', 'fasta']:
        input_fasta = input_file
    else:
        input_fasta = '%s/pos.fa' % out_dir
        cmd = 'bedtools getfasta -fi %s/assembly/%s.fa -bed %s -fo %s' % (os.environ[genome], genome.lower(), input_file, input_fasta)
        subprocess.call(cmd, shell=True)
    return input_fasta


def shuffle_input(input_file, genome, out_dir):
    ''' Shuffle the input BED and return a FASTA file '''

    if genome == 'HG19':
        genome_fasta = '%s/assembly/hg19.fa' % os.environ['HG19']
        genome_file = '%s/assembly/human.hg19.genome' % os.environ['HG19']
        gaps_file = '%s/assembly/hg19_gaps.bed' % os.environ['HG19']
    else:
        print('%s assembly files not hard-coded for shuffle' % genome, file=sys.stderr)

    shuffle_bed = '%s/neg.bed' % out_dir
    cmd = 'bedtools shuffle -excl %s -i %s -g %s > %s' % (gaps_file, input_file, genome_file, shuffle_bed)
    subprocess.call(cmd, shell=True)

    shuffle_fasta = '%s/neg.fa' % out_dir
    cmd = 'bedtools getfasta -fi %s -bed %s -fo %s' % (genome_fasta, shuffle_bed, shuffle_fasta)
    subprocess.call(cmd, shell=True)

    return shuffle_fasta

def trim_parens(string):
    if string[0] == '(':
        string = string[1:]
    if string[-1] == ')':
        string = string[:-1]
    return string

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
