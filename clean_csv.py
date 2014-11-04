#!/usr/bin/env python
from optparse import OptionParser

################################################################################
# clean_csv.py
#
# Clean up an excel-saved .csv file with \r's and commas.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    file_str = open(args[0]).readline()
    for line in file_str.split('\r'):
        a = line.split(',')
        print '\t'.join(a)
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
