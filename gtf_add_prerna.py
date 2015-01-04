#!/usr/bin/env python
from optparse import OptionParser
import copy, os, subprocess, sys, tempfile
import gff

################################################################################
# gtf_add_prerna.py
#
# Add spans of isoforms to a GTF file
#
# WARNING: Output GTF is unsorted.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <ref_gtf> <prerna_gtf>'
    parser = OptionParser(usage)
    parser.add_option('-m', dest='max_genes_overlapped', default=None, type='int', help='Don\'t include isoforms that overlap more than this many genes [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide reference GTF and output prerna GTF')
    else:
        ref_gtf = args[0]
        prerna_gtf = args[1]

    # read transcripts for filtering/processing
    transcripts = gff.read_genes(ref_gtf, key_id='transcript_id')

    # add unspliced single exon transcripts to hash
    prerna_hash = set()
    for tid in transcripts:
        tx = transcripts[tid]
        if len(tx.exons) == 1:
            tx_key = (tx.chrom, tx.exons[0].start, tx.exons[0].end, tx.strand)
            prerna_hash.add(tx_key)

    # process transcripts
    prerna_out = open(prerna_gtf, 'w')
    prerna_index = 0
    for tid in transcripts:
        tx = transcripts[tid]
        pre_start = tx.exons[0].start
        pre_end = tx.exons[-1].end
        pre_key = (tx.chrom, pre_start, pre_end, tx.strand)

        # print exons
        for i in range(len(tx.exons)):
            cols = (tx.chrom, 'dk', 'exon', str(tx.exons[i].start), str(tx.exons[i].end), '.', tx.strand, '.', gff.kv_gtf(tx.kv))
            print >> prerna_out, '\t'.join(cols)

        # print prernas
        if not pre_key in prerna_hash:
            prerna_hash.add(pre_key)
            pre_kv = copy.copy(tx.kv)
            pre_kv['transcript_id'] = 'PRERNA%d' % prerna_index
            pre_kv['transcript_type'] = 'prerna'
            prerna_index += 1
            cols = (tx.chrom, 'dk', 'exon', str(pre_start), str(pre_end), '.', tx.strand, '.', gff.kv_gtf(pre_kv))
            print >> prerna_out, '\t'.join(cols)

    prerna_out.close()

    if options.max_genes_overlapped != None:
        # intersect with self and compute overlap sets
        p = subprocess.Popen('intersectBed -wo -s -a %s -b %s' % (prerna_gtf, prerna_gtf), shell=True, stdout=subprocess.PIPE)

        tx_overlaps = {}
        for line in p.stdout:
            a = line.split('\t')

            kv1 = gff.gtf_kv(a[8])
            tid1 = kv1['transcript_id']

            if tid1.startswith('PRERNA'):
                gid1 = kv1['gene_id']
                gid2 = gff.gtf_kv(a[17])['gene_id']

                if gid1 != gid2:
                    tx_overlaps.setdefault(tid1,set()).add(gid2)

        p.communicate()

        # filter into a temp gtf
        prerna_tmp_fd, prerna_tmp_file = tempfile.mkstemp()
        prerna_out = open(prerna_tmp_file, 'w')
        for line in open(prerna_gtf):
            a = line.split('\t')
            kv = gff.gtf_kv(a[8])
            tid = kv['transcript_id']
            if len(tx_overlaps.get(tid,[])) <= options.max_genes_overlapped:
                print >> prerna_out, line,
        prerna_out.close()

        # rewrite temp to the final output
        prerna_out = open(prerna_gtf, 'w')
        for line in open(prerna_tmp_file):
            print >> prerna_out, line,
        prerna_out.close()

        os.close(prerna_tmp_fd)
        os.remove(prerna_tmp_file)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
