#!/usr/bin/env python
from optparse import OptionParser
import gff

################################################################################
# gtf_span.py
#
# Merge all of the transcripts in one gene into a single spanning gtf entry.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf_file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    gtf_file = args[0]

    genes = {}

    for line in open(gtf_file):
        a = line.split()
        gene_id = a[9][1:-2]
        genes.setdefault(gene_id,[]).append(line)

    for gene_id in genes:
        start = min([int(line.split()[3]) for line in genes[gene_id]])
        end = max([int(line.split()[4]) for line in genes[gene_id]])
        
        a = genes[gene_id][0].split('\t')
        kv = gff.gtf_kv(a[8])
        succinct_kv = {'gene_id':kv['gene_id']}
        succinct_kv['transcript_id'] = ','.join(list(set([line.split()[11][1:-2] for line in genes[gene_id]])))
        
        d = [a[0], 'gtf', 'gene', str(start), str(end), '.', a[6], '.', gff.kv_gtf(succinct_kv)]
        print '\t'.join(d)

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
