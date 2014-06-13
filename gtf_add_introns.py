#!/usr/bin/env python
from optparse import OptionParser
import copy, os, subprocess, tempfile

################################################################################
# gtf_add_introns.py
#
# Add exon-intron-exon isoforms to a GTF file, and return only exons.
#
# WARNING: Output GTF is unsorted.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <ref_gtf>'
    parser = OptionParser(usage)
    parser.add_option('-e', dest='exons_adjacent', default=False, action='store_true', help='Include adjacent exons with every intron isoform [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide GTF file')
    else:
        ref_gtf = args[0]

    transcripts = read_genes(ref_gtf, key_id='transcript_id')
    g2t_map = g2t(ref_gtf)

    ############################################
    # make new gtf
    ############################################
    intron_index = 0

    for gene_id in g2t_map:
        # make introns
        raw_introns = set()
        for transcript_id in g2t_map[gene_id]:
            if transcript_id in transcripts:
                tx = transcripts[transcript_id]
                for i in range(len(tx.exons)-1):
                    if options.exons_adjacent:
                        istart = tx.exons[i].start
                        iend = tx.exons[i+1].end
                    else:
                        istart = tx.exons[i].end+1
                        iend = tx.exons[i+1].start-1

                    raw_introns.add((istart,iend))

        # filter highly redundant intron isoforms
        introns = filter_introns(raw_introns)

        # print exon isoforms
        for transcript_id in g2t_map[gene_id]:
            if transcript_id in transcripts:
                tx = transcripts[transcript_id]
                for i in range(len(tx.exons)):
                    cols = (tx.chrom, 'dk', 'exon', str(tx.exons[i].start), str(tx.exons[i].end), '.', tx.strand, '.', kv_gtf(tx.kv))
                    print '\t'.join(cols)

        # print intron isoforms
        for istart, iend in introns:
            pre_kv = copy.copy(tx.kv)
            pre_kv['transcript_id'] = 'INTRON%d' % intron_index
            pre_kv['transcript_type'] = 'intron'
            intron_index += 1
            cols = (tx.chrom, 'dk', 'exon', str(istart), str(iend), '.', tx.strand, '.', kv_gtf(pre_kv))
            print '\t'.join(cols)


################################################################################
# filter_introns
#
# Collapse clusters of highly similar introns.
################################################################################
def filter_introns(raw_introns, overlap_diff=5):
    # initialize a cluster for every intron
    cluster_map = {}
    for i in range(len(raw_introns)):
        cluster_map[i] = i

    clusters = {}
    for i in range(len(raw_introns)):
        clusters[i] = set([i])

    # recklessly merge intron clusters
    introns_list = list(raw_introns)
    for i in range(len(raw_introns)):
        for j in range(i+1,len(raw_introns)):
            start_i, end_i = introns_list[i]
            start_j, end_j = introns_list[j]

            overlap = min(end_i, end_j) - max(start_i, start_j) + 1
            if overlap > 0:
                diff_i = abs(overlap - (end_i - start_i + 1))
                diff_j = abs(overlap - (end_j - start_j + 1))

                if diff_i <= overlap_diff and diff_j <= overlap_diff:
                    cluster_j = cluster_map[j]

                    # move cluster j to cluster i
                    for k in clusters[cluster_j]:
                        # add it to i's cluster
                        clusters[cluster_map[i]].add(k)
                        # re-map it to i's cluster
                        cluster_map[k] = cluster_map[i]

                    # erase cluster j
                    clusters[cluster_j] = set()

    # return one intron per cluster
    final_introns = set()
    for i in range(len(raw_introns)):
        if len(clusters[i]) > 0:
            final_introns.add(introns_list[list(clusters[i])[0]])

    return final_introns


################################################################################
# gtf_kv
#
# Convert the last gtf section of key/value pairs into a dict.
################################################################################
def gtf_kv(s):
    d = {}

    a = s.split(';')
    for key_val in a:
        if key_val.strip():
            eq_i = key_val.find('=')
            if eq_i != -1 and key_val[eq_i-1] != '"':
                kvs = key_val.split('=')
            else:
                kvs = key_val.split()

            if len(kvs) == 2:
                key = kvs[0]
                if kvs[1][0] == '"' and kvs[1][-1] == '"':
                    val = kvs[1].strip()[1:-1]
                else:
                    val = kvs[1].strip()
                d[key] = val

    return d


################################################################################
# g2t
#
# Given a gtf file, return a mapping of gene_id's to sets of transcript_id's
################################################################################
def g2t(gtf_file):
    d = {}

    gtf_in = open(gtf_file)

    # ignore header
    line = gtf_in.readline()
    while line[:2] == '##':
        line = gtf_in.readline()

    for line in gtf_in:
        a = line.split('\t')
        kv = gtf_kv(a[8])
        d.setdefault(kv['gene_id'],set()).add(kv['transcript_id'])

    return d


################################################################################
# kv_gtf
#
# Convert a kv hash to str gtf representation.
################################################################################
def kv_gtf(d):
    s = ''
    for key in sorted(d.keys()):
        s += '%s "%s"; ' % (key,d[key])
    return s


################################################################################
# read_genes
#
# Parse a gtf file and return a set of Gene objects in a hash keyed by the
# id given.
################################################################################
def read_genes(gtf_file, key_id='transcript_id', sort=True):
    genes = {}

    gtf_in = open(gtf_file)

    # ignore header
    line = gtf_in.readline()
    while line[:2] == '##':
        line = gtf_in.readline()
    
    while line:
        a = line.split('\t')
        if a[2] == 'exon':
            kv = gtf_kv(a[8])
            if not kv[key_id] in genes:
                genes[kv[key_id]] = Gene(a[0], a[6], kv)

            genes[kv[key_id]].add_exon(int(a[3]), int(a[4]), sort=sort)

        line = gtf_in.readline()

    gtf_in.close()

    return genes


################################################################################
# Gene
################################################################################
class Gene:
    def __init__(self, chrom, strand, kv):
        self.chrom = chrom
        self.strand = strand
        self.kv = kv
        self.exons = []

    def add_exon(self, start, end, sort=True):
        self.exons.append(Exon(start,end))
        if sort and len(self.exons) > 1 and self.exons[-2].end > start:
            #print >> sys.stderr, 'Warning: exons are not sorted - %s' % kv_gtf(self.kv)
            self.exons.sort()

    def __str__(self):
        return '%s %s %s %s' % (self.chrom, self.strand, kv_gtf(self.kv), ','.join([ex.__str__() for ex in self.exons]))

################################################################################
# Exon
################################################################################
class Exon:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __cmp__(self, x):
        if self.start < x.start:
            return -1
        elif self.start > x.start:
            return 1
        else:
            return 0

    def __str__(self):
        return 'exon(%d-%d)' % (self.start,self.end)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
