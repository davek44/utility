#!/usr/bin/env python
from optparse import OptionParser
import os, glob

'''
clear_slurm

Helper script to clear out slurm log files.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <file_size>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    file_size = None
    if len(args) > 0:
        file_size = int(args[0])

    for slurm_out in glob.glob('slurm*.out'):
        if file_size is None or os.path.getsize(slurm_out) == file_size:
            os.remove(slurm_out)
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
