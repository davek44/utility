#!/usr/bin/env python
from optparse import OptionParser
import random
import sys

################################################################################
# reservoir_sample.py
#
# Randomly choose a subset of lines in a file using single pass
# reservoir sampling.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <sample_num> <input_file>'
    parser = OptionParser(usage)
    parser.add_option('-d', dest='header', default=False, action='store_true')
    (options,args) = parser.parse_args()

    if len(args) != 2:
    	parser.error('Must provide file and sample number')
    else:
    	sample_num = int(args[0])
    	input_file = args[1]

    reservoir = ['']*sample_num

    if input_file in ['-','stdin']:
        input_in = sys.stdin
    else:
        input_in = open(input_file)

    if options.header:
        print(input_in.readline(), end='')

    # fill
    i = 0
    while i < sample_num:
    	reservoir[i] = input_in.readline()
    	i += 1

    # sample
    for line in input_in:
    	j = random.randint(0, i+1)
    	if j < sample_num:
    		reservoir[j] = line
    	i += 1

    # print
    print(''.join(reservoir), end='')

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
