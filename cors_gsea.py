#!/usr/bin/env python
from optparse import OptionParser
from scipy.stats import norm
import math, networkx, pdb, random
import util, fdr

################################################################################
# cors_gsea.py
#
# Run GSEA on gene expression correlations between a single gene to be annotated
# and all other genes.
#
# The main input file of correlations is just gene names and correlations one
# per line.
#
# This version forms a null distribution by examining the correlations from
# a randomly selected set of genes for each gene set.
#
# Note, however that I also don't love this method. I should be more careful
# about the randomly chosen genes being in the gene sets. Also, lincRNA
# expression patterns tend to differ from protein expression patterns. Finally,
# I think I should do pearson correlation of the logs rather than spearman
# correlation.
################################################################################

go_dir = '/Users/dk/research/common/data/gene_ontology'

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <cors file> <null file>'
    parser = OptionParser(usage)
    parser.add_option('--go_min', dest='go_min', type='int', default=10, help='Minimum number of genes assigned a GO term to consider enrichment of that term')
    parser.add_option('--go_max', dest='go_max', type='int', default=300, help='Maximum number of genes assigned a GO term to consider enrichment of that term')
    parser.add_option('-n', dest='null_samples', type='int', default=50, help='Number of null samples to obtain p-value [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide correlations file')
    else:
        cors_file = args[0]
        null_file = args[1]

    # get genes, correlations
    correlations_genes = []
    genes = []
    for line in open(cors_file):
        a = line.split()
        correlations_genes.append((abs(float(a[1])),a[0]))
        genes.append(a[0])
    correlations_genes.sort(reverse=True)

    # GO
    go_map, go_descs = read_go(set(genes))
    consider_go = [go_term for go_term in go_map if options.go_min <= len(go_map[go_term]) <= options.go_max]

    # compute null GO term enrichments
    null_go_enrichments = {}
    header = ''
    samples = 0
    for line in open(null_file):
        if line[0] == '>':
            if header:
                if samples < options.null_samples:
                    process_null_sample(null_correlations_genes, null_go_enrichments, go_map, consider_go)
                    samples += 1
            header = line[1:].rstrip()
            null_correlations_genes = []            
        else:
            a = line.split()
            null_correlations_genes.append((abs(float(a[1])),a[0]))
    if samples < options.null_samples:
        process_null_sample(null_correlations_genes, null_go_enrichments, go_map, consider_go)
        samples += 1

    for go_term in consider_go:
        null_out = open('null_%s_%d.txt' % (go_term,len(go_map[go_term])), 'w')
        print >> null_out, '\n'.join([str(enr) for enr in null_go_enrichments[go_term]])
        null_out.close()

    # do stats
    output_cols = []
    for go_term in consider_go:
        # compute enrichment
        enrichment = gsea_enrichment(correlations_genes, go_map[go_term])

        # compute p-value using normal approximation
        #p_val = (1+len([e for e in geneset_size_enrichments[go_size] if e >= enrichment])) / float(options.num_shuffles)
        (mean, sd) = util.mean_sd(null_go_enrichments[go_term])
        p_val = 1.0 - norm.cdf(enrichment, loc=mean, scale=sd)

        # output
        output_cols.append([go_term, enrichment, p_val, 99, len(go_map[go_term]), go_descs[go_term]])

    # FDR multiple hypothesis correction
    p_values = [oc[2] for oc in output_cols]
    q_values = fdr.storey(p_values)
    for i in range(len(output_cols)):
        output_cols[i][3] = q_values[i]

    for oc in output_cols:
        print '%-12s %.3f %.3e %.3e %4d %s' % tuple(oc)


################################################################################
# add_all_terms
#
# Add the term and its parents up the GO dag to the map.
################################################################################
def add_all_terms(go_map, go_graph, gene, source_term):
    #go_map.setdefault(gene,set()).add(source_term)
    go_map.setdefault(source_term,set()).add(gene)
    for term_edge in networkx.algorithms.traversal.breadth_first_search.bfs_edges(go_graph, source_term):
        #go_map[gene].add(term_edge[1])
        go_map.setdefault(term_edge[1],set()).add(gene)


################################################################################
# genes_from_cdf
#
# Get the gene names from the cdf file
################################################################################
def genes_from_cdf():
    genes = []
    for line in open('genes_lnc.cdf'):
        if line.startswith('Name=') and line[5:9] != 'NONE' and line[5:12] != 'HG-U133':
            genes.append(line[5:].rstrip()[:-3])
    return genes


################################################################################
# gsea_enrichment
#
# Compute enrichment score as defined in Subramanian et al 2005.
################################################################################
def gsea_enrichment(correlations_genes, term_genes):
    p_hit = [0.0]*(1+len(correlations_genes))
    p_miss = [0.0]*(1+len(correlations_genes))

    i = 1
    for (cor,gene) in correlations_genes:
        if gene in term_genes:
            p_hit[i] = p_hit[i-1] + cor
            p_miss[i] = p_miss[i-1]
        else:
            p_hit[i] = p_hit[i-1]
            p_miss[i] = p_miss[i-1] + 1.0
        i += 1

    p_hit = [p/p_hit[-1] for p in p_hit]
    p_miss = [p/p_miss[-1] for p in p_miss]

    return max([p_hit[i]-p_miss[i] for i in range(len(p_hit))])


################################################################################
# interpolate_null
#
# Compute a mean and sd for the null distribution by interpolating between
# the distributions evenly spaced among gene set sizes.
################################################################################
def interpolate_null(geneset_size_enrichments, size_skip, go_size):
    # find range
    size_less = go_size
    size_more = go_size
    for d in range(size_skip+1):
        if not size_less in geneset_size_enrichments:
            size_less -= 1
        if not size_more in geneset_size_enrichments:
            size_more += 1

    # compute interpolation weights
    max_dist = 1+max(go_size-size_less, size_more-go_size)
    w_less = float(max_dist - (go_size-size_less))
    w_more = float(max_dist - (size_more-go_size))

    # compute mean, sd
    mean = (w_less*util.mean(geneset_size_enrichments[size_less]) + w_more*util.mean(geneset_size_enrichments[size_more])) / (w_less+w_more)
    sd = (w_less*util.sd(geneset_size_enrichments[size_less]) + w_more*util.sd(geneset_size_enrichments[size_more])) / (w_less+w_more)

    return mean, sd


################################################################################
# load_go_graph
#
# Build a reversed directed graph from the GO dag for easy parent searching.
################################################################################
def load_go_graph():
    go_descs = {}
    go_graph = networkx.DiGraph()

    for line in open('%s/gene_ontology_ext.obo.obo' % go_dir):
        a = line.split()
        if a:
            if a[0] == 'id:':
                term_id = a[1]

            elif a[0] == 'name:':
                go_descs[term_id] = line[6:].rstrip()

            elif a[0] == 'is_a:':
                parent = a[1]
                go_graph.add_edge(term_id, parent)

            elif a[0] == 'is_obsolete:':
                del go_descs[term_id]

    return go_graph, go_descs


################################################################################
# make_null_dist
#
# Make the null distributions for various gene set sizes, spacing them out
################################################################################
def make_null_dist(go_min, go_max, num_shuffles, size_skip, go_map, correlations_genes, genes):
    # find set sizes
    set_sizes = set()
    for go_term in go_map:
        go_size = len(go_map[go_term])
        if go_min <= go_size <= go_max:
            set_sizes.add(go_size)

    # distribute null distributions evenly in set range
    geneset_size_enrichments = {}
    ss_max = max(set_sizes)
    ss = min(set_sizes)
    while ss < ss_max:
        ss_use = False
        for x in range(ss-size_skip, ss+size_skip+1):
            if x in set_sizes:
                ss_use = True
        if ss_use:
            geneset_size_enrichments[ss] = [gsea_enrichment(correlations_genes, random.sample(genes,ss)) for i in range(num_shuffles)]
        ss += size_skip

    # get the max
    if not ss_max in geneset_size_enrichments:
        geneset_size_enrichments[ss_max] = [gsea_enrichment(correlations_genes, random.sample(genes,ss_max)) for i in range(num_shuffles)]

    return geneset_size_enrichments


################################################################################
# process_null_sample
#
# Compute GO term enrichments for the given null sample of correlations and
# add them to the data structure.
################################################################################
def process_null_sample(correlations_genes, go_enrichments, go_map, consider_go):
    correlations_genes.sort(reverse=True)
    for go_term in consider_go:
        enrichment = gsea_enrichment(correlations_genes, go_map[go_term])
        go_enrichments.setdefault(go_term,[]).append(enrichment)


################################################################################
# read_go
#
# Read the mapping between RefSeq id's and GO terms, restricted to protein-
# coding genes that we have values for.
#
# go_map: GO term -> set(genes)
# go_descs: GO term -> string description
################################################################################
def read_go(genes):
    go_map = {}

    go_graph, go_descs = load_go_graph()

    for line in open('%s/refseq2go.txt' % go_dir):
        a = line.split()

        if a[0][:2] != 'NM':
            continue

        if a[3] != '-':
            continue

        gene = a[0][:a[0].find('.')]

        if gene in genes:
            add_all_terms(go_map, go_graph, gene, a[1])
            #go_descs[a[1]] = ' '.join(a[4:-2])

    return go_map, go_descs


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
