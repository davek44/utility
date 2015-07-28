#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess

################################################################################
# rm_nonxs.py
#
# Check for BAM files where I changed the XS tags, but left the original, 
# and delete the original.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    xs_bams_str = subprocess.check_output('find . -name accepted_hits_xs.bam', shell=True).strip()
    xs_bams = xs_bams_str.split('\n')

    for xs_bam in xs_bams:
        bam = xs_bam.replace('_xs','')
        if os.path.isfile(bam):
            print 'rm %s' % bam
            os.remove(bam)
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
