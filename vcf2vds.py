#!/usr/bin/env python
from optparse import OptionParser
import os
import shutil
from hail import *

'''
Name

Description...
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <vcf_file> <vds_file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide VCF and VDS files')
    else:
        vcf_file = args[0]
        vds_file = args[1]

    if os.path.isdir(vds_file):
        shutil.rmtree(vds_file)

    hc = HailContext()
    hc.import_vcf(vcf_file).write(vds_file)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
