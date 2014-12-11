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

    file_in = open(args[0])

    file_str = file_in.readline()

    if file_str.find('\r') != -1:
        for line in file_str.split('\r'):
            a = line.split(',')
            print '\t'.join(a)

    else:
        line = file_str
        while line:
            a = line.split(',')
            a[-1] = a[-1].rstrip()
            print '\t'.join(a)
            line = file_in.readline()

    file_in.close()
        
    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
