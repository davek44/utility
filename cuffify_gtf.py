#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess, tempfile
import gff

################################################################################
# cuffify_gtf.py
#
# Process the raw GENCODE GTF for Cufflinks use. This includes the following:
#  (1) Remove small RNAs.
#  (2) Remove small isoforms.
#  (3) Add tss_id
#  (4) Add p_id
#  (5) Keep only exon lines.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gencode_gtf>'
    parser = OptionParser(usage)
    parser.add_option('-l', dest='min_transcript_length', default=50, type='int')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide GENCODE GTF')
    else:
        full_gtf = args[0]

    ############################################################
    # remove small rna (and non-exon)
    ############################################################
    small_rnas = set(['miRNA','misc_RNA','snRNA','snoRNA','rRNA','Mt_rRNA'])
    sansrna_gtf_fd, sansrna_gtf_file = tempfile.mkstemp()
    sansrna_gtf_out = open(sansrna_gtf_file, 'w')

    # ignore header
    full_gtf_in = open(full_gtf)
    line = full_gtf_in.readline()
    while line[:2] == '##':
        line = full_gtf_in.readline()

    while line:
        a = line.split('\t')

        if a[2] == 'exon':
            kv = gff.gtf_kv(a[8])
            if kv['transcript_type'] not in small_rnas:
                print >> sansrna_gtf_out, line,

        line = full_gtf_in.readline()

    sansrna_gtf_out.close()

    ############################################################
    # remove tiny (unestimatable) transcripts
    ############################################################
    transcript_lengths = {}
    for line in open(sansrna_gtf_file):
        a = line.split('\t')
        if a[2] == 'exon':
            transcript_id = gff.gtf_kv(a[8])['transcript_id']
            transcript_lengths[transcript_id] = transcript_lengths.get(transcript_id,0) + int(a[4])-int(a[3])+1

    sanstiny_gtf_fd, sanstiny_gtf_file = tempfile.mkstemp()
    sanstiny_gtf_out = open(sanstiny_gtf_file, 'w')
    
    for line in open(sansrna_gtf_file):
        a = line.split('\t')
        kv = gff.gtf_kv(a[8])
        if transcript_lengths[kv['transcript_id']] >= options.min_transcript_length:
            print >> sanstiny_gtf_out, line,

    sanstiny_gtf_out.close()

    ############################################################
    # run cuffcompare to get id's
    ############################################################
    subprocess.call('cuffcompare -s $HG19/sequence/hg19.fa -CG -r %s %s' % (sanstiny_gtf_file, sanstiny_gtf_file), shell=True)

    # hash id's by oId
    tss_id = {}
    p_id = {}
    for line in open('cuffcmp.combined.gtf'):
        a = line.split('\t')
        kv = gff.gtf_kv(a[8])

        tss_id[kv['oId']] = kv['tss_id']
        if 'p_id' in kv:
            p_id[kv['oId']] = kv['p_id']

    ############################################################
    # add id's and print
    ############################################################
    unsorted_gtf_fd, unsorted_gtf_file = tempfile.mkstemp()
    unsorted_gtf_out = open(unsorted_gtf_file, 'w')

    for line in open(sanstiny_gtf_file):
        a = line.split('\t')
        kv = gff.gtf_kv(a[8])
        
        kv['tss_id'] = tss_id[kv['transcript_id']]
        if kv['transcript_id'] in p_id:
            kv['p_id'] = p_id[kv['transcript_id']]

        a[8] = gff.kv_gtf(kv)
        print >> unsorted_gtf_out, '\t'.join(a)

    unsorted_gtf_out.close()

    ############################################################
    # might as well sort it!
    ############################################################
    subprocess.call('sortBed -i %s' % unsorted_gtf_file, shell=True)

    ############################################################
    # clean
    ############################################################
    # temp
    os.close(sansrna_gtf_fd)
    os.remove(sansrna_gtf_file)
    os.close(sanstiny_gtf_fd)
    os.remove(sanstiny_gtf_file)
    os.close(unsorted_gtf_fd)
    os.remove(unsorted_gtf_file)

    # cuffcompare
    os.remove('cuffcmp.tracking')
    os.remove('cuffcmp.loci')
    os.remove('cuffcmp.combined.gtf')
    os.remove('cuffcmp.stats')


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
