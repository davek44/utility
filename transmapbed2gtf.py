#!/usr/bin/env python
from optparse import OptionParser
import gff

################################################################################
# transmapbed2gtf.py
#
# Convert the bed file that you get from the TransMap pipeline to a gtf file
# where adjacent blocks are merged.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bed file>'
    parser = OptionParser(usage)
    parser.add_option('-g', dest='orig_gtf', help='The original gtf file of the TransMap\'d genes to be used to transfer gene id\'s')
    parser.add_option('-m', dest='merge_dist', type='int', default=30, help='Minimum distance two exons can be apart for them to be merged [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide bed file')
    else:
        bed_file = args[0]

    # map transcript id's to gene id's if possible
    t2g = {}
    if options.orig_gtf:
        for line in open(options.orig_gtf):
            a = line.split('\t')
            kv = gff.gtf_kv(a[8])
            t2g[kv['transcript_id']] = kv['gene_id']

    # hash to disambiguate multi-mapping transcripts
    transcript_maps = {}

    for line in open(bed_file):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        tid = a[3]
        gid = t2g.get(a[3],a[3])

        transcript_maps[tid] = transcript_maps.get(tid,0) + 1
        if transcript_maps[tid] > 1:
            gid += '_v%d' % transcript_maps[tid]
            tid += '_v%d' % transcript_maps[tid]

        gene_start = int(a[1])
        gene_end = int(a[2])

        block_sizes = [int(x) for x in a[10].split(',') if x]
        block_starts = [int(x) for x in a[11].split(',') if x]

        exon_cols = []
        last_end = None
        exon_num = 1
        for i in range(len(block_starts)):
            exon_start = gene_start+1+block_starts[i]
            exon_end = gene_start+1+block_starts[i]+block_sizes[i]-1

            if last_end and last_end+options.merge_dist >= exon_start:
                # merge w/ last
                exon_cols[-1][4] = str(exon_end)
            else:
                exon_cols.append([a[0], 'TransMap', 'exon', str(exon_start), str(exon_end), '.', a[5], '.', 'gene_id "%s"; transcript_id "%s"; exon_number "%d"' % (gid,tid,exon_num)])
                exon_num += 1
            
            last_end = exon_end

        for cols in exon_cols:
            print '\t'.join(cols)
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
