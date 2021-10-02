#!/usr/bin/env python
from optparse import OptionParser
from collections import OrderedDict
import os
import pdb
import shutil
import sys

import numpy as np
from scipy.sparse import csc_matrix, dok_matrix
import pysam

import util

'''
fimo_vcf.py

Score a VCF file for motif changes.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <fasta_file> <vcf_file>'
    parser = OptionParser(usage)
    parser.add_option('-d', dest='database',
    		default='%s/jaspar/jaspar2020_core_meme.txt' % os.environ['HG38'])
    parser.add_option('-m', dest='motif',
    		default=None, help='Motif subset [Default: %default]')
    parser.add_option('-o', dest='out_dir',
    		default='fimo_vcf')
    parser.add_option('-p', dest='pvalue_t',
    		default=1e-4, type='float')
    parser.add_option('-w', dest='width',
    		default=25, type='int',
    		help='Width to search for motifs [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
    	parser.error('')
    else:
    	fasta_file = args[0]
    	vcf_file = args[1]

    if not os.path.isdir(options.out_dir):
    	os.mkdir(options.out_dir)

    assert(options.width % 2 == 1)
    half_width = options.width // 2

    ##################################################################
    # allele fasta

    fasta_open = pysam.Fastafile(fasta_file)

    fasta_ref_file = '%s/ref.fa' % options.out_dir
    fasta_ref_out = open(fasta_ref_file, 'w')
    fasta_alt_file = '%s/alt.fa' % options.out_dir
    fasta_alt_out = open(fasta_alt_file, 'w')

    snp_indexes = OrderedDict()
    si = 0

    for line in open(vcf_file):
    	if line[0] != '#':
    		a = line.split('\t')
    		chrm = a[0]
    		pos = int(a[1])
    		snp_id = a[2]
    		ref_nt = a[3]
    		alt_nt = a[4]
    		if len(ref_nt) > 1 or len(alt_nt) > 1:
    			print('Indels not implemented.', file=sys.stderr)
    			exit(1)
    		else:
    			seq_start = pos-1-half_width
    			seq_end = seq_start+options.width
    			seq_ref = fasta_open.fetch(chrm, seq_start, seq_end)

    			snp_indexes[snp_id] = si
    			si += 1

    			seq_ref_nt = seq_ref[half_width]
    			if seq_ref_nt != ref_nt:
    				print('FASTA ref %s does not match VCF ref %s' % (seq_ref_nt,ref_nt), file=sys.stderr)
    				exit(1)
    			else:
    				seq_alt = seq_ref[:half_width] + alt_nt + seq_ref[half_width+1:]
    				assert(len(seq_ref) == len(seq_alt))

    				print('>%s\n%s' % (snp_id,seq_ref), file=fasta_ref_out)
    				print('>%s\n%s' % (snp_id,seq_alt), file=fasta_alt_out)

    fasta_ref_out.close()
    fasta_alt_out.close()

    fasta_open.close()

    ##################################################################
    # FIMO

    # clean dir
    fimo_ref_dir = '%s/fimo_ref' % options.out_dir
    if os.path.isdir(fimo_ref_dir):
    	shutil.rmtree(fimo_ref_dir)
    fimo_alt_dir = '%s/fimo_alt' % options.out_dir
    if os.path.isdir(fimo_alt_dir):
    	shutil.rmtree(fimo_alt_dir)

    fimo_pvalue_t = 10*options.pvalue_t
    fimo_opts = '--thresh %e' % fimo_pvalue_t
    if options.motif is not None:
    	fimo_opts += ' --motif %s' % options.motif

    cmd_ref = 'fimo %s -o %s %s %s 2> %s.err' % (fimo_opts, fimo_ref_dir, options.database, fasta_ref_file, fimo_ref_dir)
    cmd_alt = 'fimo %s -o %s %s %s 2> %s.err' % (fimo_opts, fimo_alt_dir, options.database, fasta_alt_file, fimo_alt_dir)
    util.exec_par([cmd_ref,cmd_alt])

    # index motifs
    if options.motif is not None:
    	motif_indexes = {options.motif: 0}
    else:
	    motif_indexes = OrderedDict()
	    mi = 0
	    for line in open(options.database):
	    	if line.startswith('MOTIF'):
	    		a = line.split()
	    		motif_indexes[a[1]] = mi
	    		mi += 1


    # read output
    ref_fimo_file = '%s/fimo.tsv' % fimo_ref_dir
    alt_fimo_file = '%s/fimo.tsv' % fimo_alt_dir
    ref_motif_score, ref_motif_nlp = read_fimo_output(ref_fimo_file, snp_indexes, motif_indexes)
    alt_motif_score, alt_motif_nlp = read_fimo_output(alt_fimo_file, snp_indexes, motif_indexes)


    ##################################################################
    # compute scores

    # convert to csc
    ref_motif_score = ref_motif_score.tocsc()
    alt_motif_score = alt_motif_score.tocsc()
    
    snp_motif_score = dok_matrix(ref_motif_score.shape, dtype='float32')
    # snp_motif_nlp = csc_matrix(ref_motif_score.shape, dtype='float32')

    nlp_t = -np.log10(options.pvalue_t)
    snp_motif_nlp = ref_motif_nlp.maximum(alt_motif_nlp)
    snp_motif_mask = (snp_motif_nlp >= nlp_t)

    for motif_id, mi in motif_indexes.items():
    	# extract motif score vectors
    	ref_scores = ref_motif_score[:,mi].toarray()
    	alt_scores = alt_motif_score[:,mi].toarray()

    	# clip to min score
    	min_score = min(ref_scores.min(), alt_scores.min())
    	ref_scores = np.clip(ref_scores, min_score, np.inf)
    	alt_scores = np.clip(alt_scores, min_score, np.inf)

    	# compute differences
    	snp_scores = alt_scores - ref_scores

    	# save p-value
    	snp_mask = snp_motif_mask[:,mi]
    	snp_nlp = snp_motif_nlp[:,mi]
    	snp_nlp_mask = snp_mask.multiply(snp_nlp)
    	snp_motif_nlp[:,mi] = snp_nlp_mask

    	# save score
    	snp_scores_mask = snp_mask.multiply(snp_scores)
    	#snp_scores_mask = np.expand_dims(snp_scores_mask,-1)
    	snp_motif_score[:,mi] = snp_scores_mask


    ##################################################################
    # output

    table_out = open('%s/table.txt' % options.out_dir, 'w')
    for motif_id, mi in motif_indexes.items():
    	snp_score = snp_motif_score[:,mi].toarray().squeeze()
    	snp_nlp = snp_motif_nlp[:,mi].toarray().squeeze()
    	snp_value = np.power(10,-snp_nlp)
    	for snp_id, si in snp_indexes.items():
    		cols = [snp_id, motif_id, '%.3f'%snp_score[si], '%.1e'%snp_value[si]]
    		print('\t'.join(cols), file=table_out)
    table_out.close()



def read_fimo_output(fimo_out_file, snp_indexes, motif_indexes):
	num_snps = len(snp_indexes)
	num_motifs = len(motif_indexes)

	snp_motif_score = dok_matrix((num_snps,num_motifs), dtype='float32')
	snp_motif_nlp = dok_matrix((num_snps,num_motifs), dtype='float32')

	fimo_out_open = open(fimo_out_file)
	fimo_out_open.readline()

	for line in fimo_out_open:
		a = line.split()
		if len(a) > 0 and line[0] != '#':
			motif_id = a[0]
			snp_id = a[2]
			score = float(a[6])
			pval = float(a[7])

			si = snp_indexes[snp_id]
			mi = motif_indexes[motif_id]

			snp_motif_score[si,mi] = max(score, snp_motif_score[si,mi])
			snp_motif_nlp[si,mi] = max(-np.log10(pval), snp_motif_nlp[si,mi])

	fimo_out_open.close()

	return snp_motif_score, snp_motif_nlp

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
