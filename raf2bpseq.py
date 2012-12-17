#!/usr/bin/env python
from optparse import OptionParser
import random, sys

################################################################################
# raf2bpseq.py
#
# Convert the RAF output to a bpseq file for visualization.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <raf file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide RAF output file')
    else:
        raf_out = args[0]

    # get multiple sequence alignment and predicted structure
    msa, structure = parse_raf(raf_out)

    # compute consensus sequence
    cons = compute_consensus(msa, structure)

    # remove gaps
    cons_nogap = []
    structure_nogap = []
    for i in range(len(cons)):
        if cons[i] != '-':
            cons_nogap.append(cons[i])
            structure_nogap.append(structure[i])

    # map nucleotides to their pairs
    nt2nt = [0]*len(cons_nogap)
    for (i,j) in get_pairs(structure_nogap):
        nt2nt[i] = str(j+1)
        nt2nt[j] = str(i+1)

    # print non-gapped sequence and structure in bpseq format
    for i in range(len(nt2nt)):
        print (i+1), cons_nogap[i], nt2nt[i]


################################################################################
# compute_consensus
#
# Compute consensus sequence using the multiple sequence alignment. Choose a
# random nt when there is a tie. Rescue gapped consensus's that are placed
# in the structure (stupidly by RAF).
################################################################################    
def compute_consensus(msa, structure):
    cons = []
    for j in range(len(msa[0])):
        nts = {}
        for i in range(len(msa)):
            nts[msa[i][j]] = nts.get(msa[i][j],0) + 1

        max_count = max(nts.values())
        max_nt = random.choice([nt for nt in nts.keys() if nts[nt] == max_count])

        # rescue base pairs in gaps
        if max_nt == '-' and structure[j] != '.':
            print >> sys.stderr, 'Consensus gap in structure base pair - %d' % (j+1)
            
            del nts['-']
            max_count = max(nts.values())
            max_nt = random.choice([nt for nt in nts.keys() if nts[nt] == max_count])

        cons.append(max_nt)

    return cons


################################################################################
# get_pairs
#
# Extract a list of paired bases from the structure list
################################################################################        
def get_pairs(structure):
    pairs = []
    sq = []
    for i in range(len(structure)):
        if structure[i] == '(':
            sq.append(i)
        elif structure[i] == ')':
            pairs.append((sq.pop(),i))
    return pairs


################################################################################
# parse_raf
#
# Return the structure and multiple sequence alignment from RAF output.
################################################################################
def parse_raf(raf_out):
    msa = []
    for line in open(raf_out):
        if line[0] == '>':
            msa.append([])
        else:
            msa[-1] += list(line.rstrip())

    return msa[:-1], msa[-1]


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
