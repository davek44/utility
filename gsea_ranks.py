#!/usr/bin/env python
from optparse import OptionParser
import gsea

################################################################################
# gsea_ranks.py
#
# Print out the ranks for a given GO term.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <cors file> <GO term>'
    parser = OptionParser(usage)
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide correlations file')
    else:
        cors_file = args[0]
        go_term = args[1]

    # get genes, correlations
    correlations_genes = []
    genes = []
    for line in open(cors_file):
        a = line.split()
        correlations_genes.append((abs(float(a[1])),a[0]))
        genes.append(a[0])
    correlations_genes.sort(reverse=True)

    # GO
    go_map, go_descs = gsea.read_go(set(genes))
    go_term_map = go_map[go_term]

    # print ranks, correlations
    i = 1
    for (cor,gene) in correlations_genes:
        if gene in go_term_map:
            print i, cor
        i += 1


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
