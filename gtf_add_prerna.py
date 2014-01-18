#!/usr/bin/env python
from optparse import OptionParser
import copy, os, subprocess, tempfile
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
    usage = 'usage: %prog [options] <ref_gtf>'
    parser = OptionParser(usage)
    parser.add_option('-e', dest='exons_adjacent', default=False, action='store_true', help='Include adjacent exons with every intron isoform [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide GTF file')
    else:
        ref_gtf = args[0]

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
    prerna_index = 0
    for tid in transcripts:
        tx = transcripts[tid]
        pre_start = tx.exons[0].start
        pre_end = tx.exons[-1].end
        pre_key = (tx.chrom, pre_start, pre_end, tx.strand)

        # print exons
        for i in range(len(tx.exons)):
            cols = (tx.chrom, 'dk', 'exon', str(tx.exons[i].start), str(tx.exons[i].end), '.', tx.strand, '.', gff.kv_gtf(tx.kv))
            print '\t'.join(cols)

        # print prernas
        if not pre_key in prerna_hash:
            prerna_hash.add(pre_key)
            pre_kv = copy.copy(tx.kv)
            pre_kv['transcript_id'] = 'PRERNA%d' % prerna_index
            pre_kv['transcript_type'] = 'prerna'
            prerna_index += 1
            cols = (tx.chrom, 'dk', 'exon', str(pre_start), str(pre_end), '.', tx.strand, '.', gff.kv_gtf(pre_kv))
            print '\t'.join(cols)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
