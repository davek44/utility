#!/usr/bin/env python
from optparse import OptionParser

################################################################################
# gtf2utrs.py
#
# Take a gtf file with exons and CDS annotated and return a gtf of the UTRs.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gtf file>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide gtf file')
    else:
        gtf_file = args[0]
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
