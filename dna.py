#!/usr/bin/env python
# no __future__ division!
import string, sys, math, random, pdb
from collections import OrderedDict

############################################################
# dna
#
# Common methods for dealing with dna sequnces
############################################################

############################################################
# fasta2dict
#
# Read a multifasta file into a dict.  Taking the whole line
# as the key.
#
# I've found this can be quite slow for some reason, even
# for a single fasta entry.
############################################################
def fasta2dict(fasta_file):
    fasta_dict = OrderedDict()
    header = ''

    for line in open(fasta_file):
        if line[0] == '>':
            #header = line.split()[0][1:]
            header = line[1:].rstrip()
            fasta_dict[header] = ''
        else:
            fasta_dict[header] += line.rstrip()

    return fasta_dict


############################################################
# rc
#
# Reverse complement sequence
############################################################
def rc(seq):
    return seq.translate(string.maketrans("ATCGatcg","TAGCtagc"))[::-1]


############################################################
# rc_file
#
# Reverse complement sequences in a file
############################################################
def rc_file(seq_file):
    seq = ''
    for line in open(seq_file):
        if line[0] == '>':
            # print last sequence
            if seq:
                print(seq.translate(string.maketrans("ATCGatcg","TAGCtagc"))[::-1])

            # print header
            print('%s_rc' % line.rstrip())
            seq = ''
        else:
            seq += line.rstrip()

    # print last sequence
    if seq:
        print(seq.translate(string.maketrans("ATCGatcg","TAGCtagc"))[::-1])


############################################################
# count_kmers
#
# Count kmers from forward and reverse strand
############################################################
def count_kmers(k, seq, all=False):
    kmers = {}
    N = len(seq)
    rc_seq = rc(seq)
    for i in range(N-k+1):
        # forward
        kmer = seq[i:i+k]
        if kmers.has_key(kmer):
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1

        # reverse
        kmer = rc_seq[i:i+k]
        if kmers.has_key(kmer):
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1

    # remove non-ACGT kmers
    nts = {'A':1, 'C':1, 'G':1, 'T':1}
    for kmer in kmers.keys():
        for nt in kmer:
            if not nts.has_key(nt):
                del kmers[kmer]
                break

    if all:
        # add zero count kmers
        for i in range(int(math.pow(4,k))):
            kmer = int2kmer(k,i)
            if not kmers.has_key(kmer):
                kmers[kmer] = 0

    return kmers


############################################################
# count_kmers_file
#
# Count kmers from forward and reverse strand
############################################################
def count_kmers_file(k, fasta_file, all=False):
    kmers = {}
    seq = None
    for line in open(fasta_file):
        if line[0] == '>':
            if seq:
                seq_kmers = count_kmers(k, seq, all)
                for kmer in seq_kmers:
                    kmers[kmer] = kmers.get(kmer,0) + seq_kmers[kmer]
            seq = ''
        else:
            seq += line.rstrip()

    seq_kmers = count_kmers(k, seq, all)
    for kmer in seq_kmers:
        kmers[kmer] = kmers.get(kmer,0) + seq_kmers[kmer]

    return kmers


############################################################
# int2kmer
#
# Map integers to kmers
############################################################
def int2kmer(k,num):
    nts = ['A','C','G','T']
    kmer = ''
    for x in range(k):
        b = int(math.pow(4, k-1-x))
        kmer += nts[num / b]
        num =  num % b
    return kmer


############################################################
# canonical_kmers
#
# Clean up a dict of kmer counts by combining kmers with
# their reverse complements.  All counts are then divided
# by 2.  Careful about palindromes.
############################################################
def canonical_kmers(kmers, return_all=False):
    canon_kmers = {}
    for kmer in kmers:
        kmer_rc = rc(kmer)

        if kmer < kmer_rc:
            # add current
            if canon_kmers.has_key(kmer):
                canon_kmers[kmer] += kmers[kmer] / 2.0
            else:
                canon_kmers[kmer] = kmers[kmer] / 2.0

        elif kmer_rc < kmer:
            # add current
            if canon_kmers.has_key(kmer_rc):
                canon_kmers[kmer_rc] += kmers[kmer] / 2.0
            else:
                canon_kmers[kmer_rc] = kmers[kmer] / 2.0

        elif kmer == kmer_rc:
            # add once, divide by 2 bc we double counted it
            #  once on each strand
            canon_kmers[kmer] = kmers[kmer] / 2.0

    if return_all:
        # add back reverse complements
        for kmer in kmers:
            if not canon_kmers.has_key(kmer):
                canon_kmers[kmer] = canon_kmers[rc(kmer)]

    return canon_kmers


############################################################
# fasta_rand
#
# Randomly sample 'num_seq' sequences from a multi-fasta
# file, with an option to draw pairs of mates.
############################################################
def fasta_rand(num_seq, reads_file, out_file='', mates_file=''):
    random.seed()

    seqs = fasta2dict(reads_file)

    if out_file:
        out = open(out_file, 'w')
    else:
        out = sys.stdout

    if mates_file:
        # get mates
        mates = {}
        for line in open(mates_file):
            (lr,rr) = line.split()
            mates[lr] = rr
            mates[rr] = lr

        # sample from left reads, print both
        for h in random.sample(mates.keys(), num_seq/2):
            print >> out, '>%s' % h
            print >> out, seqs[h]
            print >> out, '>%s' % mates[h]
            print >> out, seqs[mates[h]]

    else:
        # sample from all
        for h in random.sample(seqs.keys(), num_seq):
            print >> out, '>%s' % h
            print >> out, seqs[h]


############################################################
# fastq_rand
#
# Randomly sample 'num_seq' sequences from a fastq file
############################################################
def fastq_rand(num_seq, reads_file, out_file=''):
    random.seed()

    # count sequences
    total_seq = 0
    fqf = open(reads_file)
    header = fqf.readline()
    while header:
        seq = fqf.readline()
        mid = fqf.readline()
        qual = fqf.readline()
        total_seq += 1
        header = fqf.readline()
    fqf.close()

    # choose random sequences
    rand_seqs = sorted(random.sample(xrange(total_seq), num_seq))

    if out_file:
        out = open(out_file, 'w')
    else:
        out = sys.stdout

    fqf = open(reads_file)
    seq_i = 0
    rand_i = 0
    header = fqf.readline()
    while header:
        seq = fqf.readline()
        mid = fqf.readline()
        qual = fqf.readline()

        if seq_i == rand_seqs[rand_i]:
            print >> out, header + seq + mid + qual,
            rand_i += 1
            if rand_i >= num_seq:
                break
        seq_i += 1
        header = fqf.readline()
    fqf.close()

############################################################
# nt_composition
#
# Return a dict of the nt counts in the given sequence,
# making no assumptions about what the sequence components
# are.
############################################################
def nt_composition(seq):
    comp = {}
    for nt in seq:
        nt = nt.upper()
        if comp.has_key(nt):
            comp[nt] += 1
        else:
            comp[nt] = 1
    return comp


############################################################
# nt_composition_file
#
# Return a dict of the nt counts of the sequences in the
# given file making no assumptions about what the sequence
# components are.
############################################################
def nt_composition_file(seq_file):
    comp = {}
    for line in open(seq_file):
        if line[0] != '>':
            for nt in line.rstrip():
                nt = nt.upper()
                if comp.has_key(nt):
                    comp[nt] += 1
                else:
                    comp[nt] = 1
    return comp

############################################################
# translate
#
# Translate a dna sequence into an amino acid.  Attempts
# to maintain lowercase or uppercase.  If a codon contains
# both lowercase and uppercase, returns a lowercase codon.
############################################################
code = {     'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C', \
             'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C', \
             'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*', \
             'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', \
             'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R', \
             'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', \
             'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', \
             'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', \
             'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', \
             'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', \
             'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', \
             'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', \
             'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', \
             'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', \
             'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', \
             'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G', \

             'ttt': 'f', 'tct': 's', 'tat': 'y', 'tgt': 'c', \
             'ttc': 'f', 'tcc': 's', 'tac': 'y', 'tgc': 'c', \
             'tta': 'l', 'tca': 's', 'taa': '*', 'tga': '*', \
             'ttg': 'l', 'tcg': 's', 'tag': '*', 'tgg': 'w', \
             'ctt': 'l', 'cct': 'p', 'cat': 'h', 'cgt': 'r', \
             'ctc': 'l', 'ccc': 'p', 'cac': 'h', 'cgc': 'r', \
             'cta': 'l', 'cca': 'p', 'caa': 'q', 'cga': 'r', \
             'ctg': 'l', 'ccg': 'p', 'cag': 'q', 'cgg': 'r', \
             'att': 'i', 'act': 't', 'aat': 'n', 'agt': 's', \
             'atc': 'i', 'acc': 't', 'aac': 'n', 'agc': 's', \
             'ata': 'i', 'aca': 't', 'aaa': 'k', 'aga': 'r', \
             'atg': 'm', 'acg': 't', 'aag': 'k', 'agg': 'r', \
             'gtt': 'v', 'gct': 'a', 'gat': 'd', 'ggt': 'g', \
             'gtc': 'v', 'gcc': 'a', 'gac': 'd', 'ggc': 'g', \
             'gta': 'v', 'gca': 'a', 'gaa': 'e', 'gga': 'g', \
             'gtg': 'v', 'gcg': 'a', 'gag': 'e', 'ggg': 'g' \
             }

def translate(dna):
    if len(dna) % 3 != 0:
        print >> sys.stderr, 'DNA sequence is not have length divisible by 3.'

    i = 0
    peptide = ''
    while i+2 < len(dna):
        if code.has_key(dna[i:i+3]):
            peptide += code[dna[i:i+3]]
        else:
            peptide += code[dna[i:i+3].lower()]
        i += 3
    return peptide

def translate_file(seq_file):
    seq = ''
    for line in open(seq_file):
        if line[0] == '>':
            # print last sequence
            if seq:
                print(translate(seq))

            # print header
            print('%s_aa' % line.rstrip())
            seq = ''
        else:
            seq += line.rstrip()

    # print last sequence
    if seq:
        print(translate(seq))

############################################################
# __main__
############################################################
if __name__ == '__main__':
    # print fasta file DNA composition
    if len(sys.argv) == 3 and sys.argv[1] == '--comp':
        nt_comp = nt_composition_file(sys.argv[2])
        nsum = float(sum(nt_comp.values()))
        at = (nt_comp['A']+nt_comp['T'])
        gc = (nt_comp['C']+nt_comp['G'])
        print('A/T %d (%.4f)' % (at,at/nsum))
        print('G/C %d (%.4f)' % (gc,gc/nsum))

    # print kmer composition
    elif len(sys.argv) == 4 and sys.argv[1] == '--kmers':
        nt_kmers = count_kmers_file(int(sys.argv[2]), sys.argv[3])
        kmer_sum = float(sum(nt_kmers.values()))
        for kmer in sorted(nt_kmers.keys()):
            print('%s %8d %.4f' % (kmer,nt_kmers[kmer],nt_kmers[kmer]/kmer_sum))

    # print DNA reverse complement
    elif len(sys.argv) == 3 and sys.argv[1] == '--rc':
        rcf = False
        for i in range(min(10,len(sys.argv[2]))):
            if 'acgtACGT'.find(sys.argv[2][i]) == -1:
                rcf = True
                break
        if rcf:
            rc_file(sys.argv[2])
        else:
            print(rc(sys.argv[2]))

    # translate DNA sequence
    elif len(sys.argv) == 3 and sys.argv[1] == '--trans':
        tf = False
        for i in range(min(10,len(sys.argv[2]))):
            if 'acgtACGT'.find(sys.argv[2][i]) == -1:
                tf = True
                break
        if tf:
            translate_file(sys.argv[2])
        else:
            print(translate(sys.argv[2]))

    # print random sequences from fasta file
    elif len(sys.argv) == 4 and sys.argv[1] == '--rand':
        fasta_rand(int(sys.argv[3]), sys.argv[2])

    # print random sequences from fastq file
    elif len(sys.argv) == 4 and sys.argv[1] == '--qrand':
        #pdb.runcall(fastq_rand, int(sys.argv[3]), sys.argv[2])
        fastq_rand(int(sys.argv[3]), sys.argv[2])

    else:
        print('Usage: dna.py --comp <seq_file>')
        print('Usage: dna.py --rc <seq>')
        print('Usage: dna.py --rand <seq_file> <num>')
        print('Usage: dna.py --kmers <k> <seq_file>')
