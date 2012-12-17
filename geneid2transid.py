#!/usr/bin/env python
from optparse import OptionParser

################################################################################
# geneid2transid.py
#
# Given a gene id, produce a transcript id to punch into the browser
################################################################################

lnc_catalog = '/Users/dk/research/common/data/lncrna/lnc_catalog.gtf'

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gene id>\nUsage: %prog [options] <gene id file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    # parse input
    if len(args) == 0:
        parser.error('Must provide gene id or file of gene ids')
    else:
        if args[0].startswith('XLOC'):
            print find_longest_transcript(args[0])
        else:
            for line in open(args[0]):
                print find_longest_transcript(line.rstrip())



################################################################################
# find_longest_transcript
#
# Return the longest transcript for this gene
################################################################################
def find_longest_transcript(gene_id):
    # find all transcripts and sum lengths
    transcripts = {}
    for line in open(lnc_catalog):
        a = line.split('\t')
        kv = gtf_kv(a[8])

        if gene_id == kv['gene_id']:
            tx = kv['transcript_id']
            transcripts[tx] = transcripts.get(tx,0) + int(a[4])-int(a[3])+1

    # return longest
    tx_len = max([l for l in transcripts.values()])
    for tx in transcripts:
        if transcripts[tx] == tx_len:
            return tx


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
            if key_val.find('=') != -1:
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
# __main__
################################################################################
if __name__ == '__main__':
    main()
