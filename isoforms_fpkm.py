#!/usr/bin/env python
from optparse import OptionParser

################################################################################
# isoforms_fpkm.py
#
# Print the FPKM values for all isoforms of the given gene.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gene_id> <iso_ft>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide a gene_id and isoforms.fpkm_tracking file')
    else:
        gene_id = args[0]
        iso_ft = args[1]

    # get headers
    fpkm_in = open(iso_ft)
    headers = fpkm_in.readline().split()

    # determine sample table length
    sample_len = 0
    for i in range(len(headers)):
        if headers[i][-5:] == '_FPKM':
            sample = headers[i][:-5]
            if len(sample) > sample_len:
                sample_len = len(sample)

    for line in fpkm_in:
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        tracking_id = a[0]
        line_gene_id = a[3]

        if line_gene_id == gene_id:
            i = 9
            while i < len(a):
                sample = headers[i][:-5]

                if a[i+3] in ['FAIL','HIDATA']:
                    cols = (tracking_id, sample_len, sample, a[i+3])
                    print '%-18s  %*s  %11s' % cols
                else:
                    fpkm = float(a[i])
                    cols = (tracking_id, sample_len, sample, fpkm)
                    print '%-18s  %*s  %11.3f' % cols

                i += 4

    fpkm_in.close()
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
