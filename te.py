#!/usr/bin/env python
from optparse import OptionParser
import glob, os, subprocess, sys
from gff import gtf_kv

################################################################################
# te.py
#
# Methods to work with transposable element annotations.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()


################################################################################
# hash_genes_repeats
#
# Hash genes in gtf_file to sets of repeats in repeats_gff.
################################################################################
def hash_genes_repeats(gtf_file, repeats_gff, gene_key='gene_id', add_star=True, stranded=False):
    gene_repeats = {}
    for line in open(gtf_file):
        a = line.split('\t')
        gene_id = gtf_kv(a[8])[gene_key]
        gene_repeats[gene_id] = set()

    p = subprocess.Popen('intersectBed -wo -a %s -b %s' % (gtf_file, repeats_gff), shell=True, stdout=subprocess.PIPE)
    line = p.stdout.readline()
    while line:
        a = line.split('\t')

        # get names
        gene_id = gtf_kv(a[8])[gene_key]
        rep_kv = gtf_kv(a[17])
        rep = rep_kv['repeat']
        fam = rep_kv['family']

        # get strands
        gene_strand = a[6]
        te_strand = a[15]

        if stranded:
            if gene_strand == te_strand:
                orient = '+'
            else:
                orient = '-'

            gene_repeats[gene_id].add((rep,fam,orient))
            if add_star:
                gene_repeats[gene_id].add(('*',fam,orient))
                gene_repeats[gene_id].add(('*','*',orient))

        else:
            gene_repeats[gene_id].add((rep,fam))
            if add_star:
                gene_repeats[gene_id].add(('*',fam))
                gene_repeats[gene_id].add(('*','*'))

        line = p.stdout.readline()
    p.communicate()

    return gene_repeats


################################################################################
# hash_genes_repeats_nt
#
# Hash genes in gtf_file to a dict mapping repeats in repeats_gff to nt overlaps.
#
# Warnings:
#  -If we hash by gene_id, we want to have chosen a single isoform per gene.
################################################################################
def hash_genes_repeats_nt(gtf_file, repeats_gff, gene_key='gene_id', add_star=True, stranded=False):
    gene_repeat_nt = {}

    p = subprocess.Popen('intersectBed -wo -a %s -b %s' % (gtf_file, repeats_gff), shell=True, stdout=subprocess.PIPE)
    line = p.stdout.readline()
    while line:
        a = line.split('\t')

        # get names
        gene_id = gtf_kv(a[8])[gene_key]
        rep_kv = gtf_kv(a[17])
        rep = rep_kv['repeat']
        fam = rep_kv['family']

        # get overlap
        nt_overlap = int(a[18])

        # get strands
        gene_strand = a[6]
        te_strand = a[15]

        if stranded:
            if gene_strand == te_strand:
                orient = '+'
            else:
                orient = '-'

            if not gene_id in gene_repeat_nt:
                gene_repeat_nt[gene_id] = {}

            gene_repeat_nt[gene_id][(rep,fam,orient)] = gene_repeat_nt[gene_id].get((rep,fam,orient),0) + nt_overlap
            if add_star:
                gene_repeat_nt[gene_id][('*',fam,orient)] = gene_repeat_nt[gene_id].get(('*',fam,orient),0) + nt_overlap
                gene_repeat_nt[gene_id][('*','*',orient)] = gene_repeat_nt[gene_id].get(('*','*',orient),0) + nt_overlap

        else:
            if not gene_id in gene_repeat_nt:
                gene_repeat_nt[gene_id] = {}

            gene_repeat_nt[gene_id][(rep,fam)] = gene_repeat_nt[gene_id].get((rep,fam),0) + nt_overlap
            if add_star:
                gene_repeat_nt[gene_id][('*',fam)] = gene_repeat_nt[gene_id].get(('*',fam),0) + nt_overlap
                gene_repeat_nt[gene_id][('*','*')] = gene_repeat_nt[gene_id].get(('*','*'),0) + nt_overlap

        line = p.stdout.readline()
    p.communicate()

    return gene_repeat_nt


################################################################################
# hash_genes_repeats_num
#
# Hash genes in gtf_file to a dict mapping repeats in repeats_gff to the number
# of repeat instances overlapped by that gene.
#
# This is tricky because if I count each overlap, I will overestimate for
# broken fragments of a single insertion. Instead, I'm normalizing overlapped
# nt by the TE length. But that won't really do the right thing for 3' biased
# L1 elements.
#
# Warnings:
#  -If we hash by gene_id, we want to have chosen a single isoform per gene.
################################################################################
def hash_genes_repeats_num(gtf_file, repeats_gff, gene_key='gene_id', add_star=True, stranded=False):
    gene_repeat_num = {}

    p = subprocess.Popen('intersectBed -wo -a %s -b %s' % (gtf_file, repeats_gff), shell=True, stdout=subprocess.PIPE)
    line = p.stdout.readline()
    while line:
        a = line.split('\t')

        # get names
        gene_id = gtf_kv(a[8])[gene_key]
        rep_kv = gtf_kv(a[17])
        rep = rep_kv['repeat']
        fam = rep_kv['family']

        # get strands
        gene_strand = a[6]
        te_strand = a[15]

        if stranded:
            if gene_strand == te_strand:
                orient = '+'
            else:
                orient = '-'

            if not gene_id in gene_repeat_num:
                gene_repeat_num[gene_id] = {}

            gene_repeat_num[gene_id][(rep,fam,orient)] = gene_repeat_num[gene_id].get((rep,fam,orient),0) + 1
            if add_star:
                gene_repeat_num[gene_id][('*',fam,orient)] = gene_repeat_num[gene_id].get(('*',fam,orient),0) + 1
                gene_repeat_num[gene_id][('*','*',orient)] = gene_repeat_num[gene_id].get(('*','*',orient),0) + 1

        else:
            if not gene_id in gene_repeat_num:
                gene_repeat_num[gene_id] = {}

            gene_repeat_num[gene_id][(rep,fam)] = gene_repeat_num[gene_id].get((rep,fam),0) + 1
            if add_star:
                gene_repeat_num[gene_id][('*',fam)] = gene_repeat_num[gene_id].get(('*',fam),0) + 1
                gene_repeat_num[gene_id][('*','*')] = gene_repeat_num[gene_id].get(('*','*'),0) + 1

        line = p.stdout.readline()
    p.communicate()

    return gene_repeat_num


################################################################################
# hash_genes_repeats_ntnum
#
# Hash genes in gtf_file to a dict mapping repeats in repeats_gff to the number
# of repeat instances overlapped by that gene.
#
# This is tricky because if I count each overlap, I will overestimate for
# broken fragments of a single insertion. Instead, I'm normalizing overlapped
# nt by the TE length. But that won't really do the right thing for 3' biased
# L1 elements.
#
# Warnings:
#  -If we hash by gene_id, we want to have chosen a single isoform per gene.
################################################################################
def hash_genes_repeats_ntnum(gtf_file, repeats_gff, gene_key='gene_id', add_star=True, stranded=False):
    # hash nt overlapped
    gene_repeat_nt = hash_genes_repeats_nt(gtf_file, repeats_gff, gene_key=gene_key, add_star=add_star, stranded=stranded)

    # determine repeat lengths
    repeat_lengths = {}
    
    # normalize
    gene_repeat_num = {}
    for gene_id in gene_repeat_nt:
        pass

    return gene_repeat_num


################################################################################
# hash_repeats_genes
#
# Hashr repeats in repeats_gff to sets of genes in gtf_file.
################################################################################
def hash_repeats_genes(gtf_file, repeats_gff, gene_key='gene_id', add_star=True, stranded=False):
    repeat_genes = {}
    if add_star:
        repeat_genes[('*','*','+')] = set()
        repeat_genes[('*','*','-')] = set()
    else:
        repeat_genes[('*','*')] = set()

    for line in open(repeats_gff):
        a = line.split('\t')
        kv = gtf_kv(a[8])
        
        if stranded:
            repeat_genes[(kv['repeat'],kv['family'],'+')] = set()
            repeat_genes[(kv['repeat'],kv['family'],'-')] = set()
            if add_star:
                repeat_genes[('*',kv['family'],'+')] = set()
                repeat_genes[('*',kv['family'],'-')] = set()
        else:
            repeat_genes[(kv['repeat'],kv['family'])] = set()
            if add_star:
                repeat_genes[('*',kv['family'])] = set()

    
    p = subprocess.Popen('intersectBed -wo -a %s -b %s' % (gtf_file, repeats_gff), shell=True, stdout=subprocess.PIPE)
    line = p.stdout.readline()
    while line:
        a = line.split('\t')

        # get names
        gene_id = gtf_kv(a[8])[gene_key]
        rep_kv = gtf_kv(a[17])
        rep = rep_kv['repeat']
        fam = rep_kv['family']

        # get strands
        gene_strand = a[6]
        te_strand = a[15]

        if stranded:
            if gene_strand == te_strand:
                orient = '+'
            else:
                orient = '-'

            repeat_genes[(rep,fam,orient)].add(gene_id)
            if add_star:
                repeat_genes[('*',fam,orient)].add(gene_id)
                repeat_genes[('*','*',orient)].add(gene_id)

        else:
            repeat_genes[(rep,fam)].add(gene_id)
            if add_star:
                repeat_genes[('*',fam)].add(gene_id)
                repeat_genes[('*','*')].add(gene_id)

        line = p.stdout.readline()
    p.communicate()

    return repeat_genes


################################################################################
# hash_repeat_family
#
# Hash repeat -> family from the RepeatMasker GFF.
################################################################################
def hash_repeat_family():
    repeat_family = {}
    for line in open('%s/hg19.fa.out.tp.gff' % os.environ['MASK']):
        a = line.split('\t')
        kv = gtf_kv(a[8])
        repeat_family[kv['repeat']] = kv['family']
    return repeat_family


################################################################################
# map_dfam_family
#
# Return a dict mapping DFAM repeats to RepeatMasker families.
################################################################################
def map_dfam_family():
    repeat_family = {}
    for line in open('%s/hg19.fa.out.tp.gff' % os.environ['MASK']):
        a = line.split('\t')
        kv = gtf_kv(a[8])
        repeat_family[kv['repeat']] = kv['family']

    dfam_family = {}
    for repeat in repeat_family:
        dfam_tes = map_rm_dfam(repeat, quiet=True)
        for dfam_te in dfam_tes:
            dfam_family[dfam_te] = repeat_family[repeat]

    return dfam_family


################################################################################
# map_dfam_repeat
#
# Return a dict mapping DFAM repeats to sets of RepeatMasker repeats.
################################################################################
def map_dfam_repeat():
    repeats = set()
    for line in open('%s/hg19.fa.out.tp.gff' % os.environ['MASK']):
        a = line.split('\t')
        kv = gtf_kv(a[8])
        repeats.add(kv['repeat'])

    dfam_repeat = {}
    for repeat in repeats:
        dfam_tes = map_rm_dfam(repeat, quiet=True)
        for dfam_te in dfam_tes:
            dfam_repeat.setdefault(dfam_te,set()).add(repeat)

    return dfam_repeat


################################################################################
# map_rm_dfam
#
# Map a RepeatMasker name to a DFAM name.
################################################################################
def map_rm_dfam(repeat, quiet=False):
    if os.path.isfile('%s/hmms/%s.hmm' % (os.environ['DFAM'],repeat)):
        dfam_reps = [repeat]
    elif os.path.isfile('%s/hmms/%sv.hmm' % (os.environ['DFAM'],repeat)):
        dfam_reps = [repeat+'v']
    else:
        hmm_files = glob.glob('%s/hmms/%s_*.hmm' % (os.environ['DFAM'],repeat))

        # if no hits
        if len(hmm_files) == 0:
            # try removing "-int"
            if repeat[-4:] == '-int' and os.path.isfile('%s/hmms/%s.hmm' % (os.environ['DFAM'],repeat[:-4])):
                dfam_reps = [repeat[:-4]]
            else:
                # missing
                if not quiet:
                    print >> sys.stderr, 'Missing DFAM name for %s' % repeat
                dfam_reps = []

        # if hits
        else:
            # grab em
            dfam_reps = []
            for i in range(len(hmm_files)):
                start = hmm_files[i].rfind('/')+1
                end = hmm_files[i].rfind('.hmm')
                dfam_reps.append(hmm_files[i][start:end])

    return dfam_reps


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
