#!/usr/bin/env python
from optparse import OptionParser
import gzip, pdb
import pysam

################################################################################
# rmdup_iclip.py
#
# Remove duplicates in Tollervey and Zamack et al's CLIP-Seq data, where the
# reads have barcodes at varying positions.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <barcode_indexes> <bam> <fastq1> ...'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) < 3:
        parser.error('Must provide barcode indexes, BAM file, and FASTQ files')
    else:
        barcode_indexes = [int(bi) for bi in args[0].split(',')]
        bam_file = args[1]
        fastq_files = args[2:]

    # map headers to barcodes
    header_barcodes = {}
    for fastq_file in fastq_files:
        if fastq_file[-2:] == 'gz':
            fastq_in = gzip.open(fastq_file)
        else:
            fastq_in = open(fastq_file)

        header = fastq_in.readline()
        while header:
            seq = fastq_in.readline()
            mid = fastq_in.readline()
            qual = fastq_in.readline()

            align_header = header[1:].split()[0]
            barcode = ''.join([seq[bi] for bi in barcode_indexes])
            header_barcodes[align_header] = barcode

            header = fastq_in.readline()

    # open BAM
    bam_in = pysam.Samfile(bam_file, 'rb')
    bam_out = pysam.Samfile(bam_file[:-4] + '_rmdup.bam', 'wb', template=bam_in)

    alignment_hash = set()

    for aligned_read in bam_in:
        # hash by chrom, start, strand, barcode
        align_key = (aligned_read.tid, aligned_read.pos, aligned_read.is_reverse, header_barcodes[aligned_read.qname])

        # if alignment not yet printed
        if not align_key in alignment_hash:
            bam_out.write(aligned_read)
            alignment_hash.add(align_key)

    bam_in.close()
    bam_out.close()

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
