#!/usr/bin/env python
from optparse import OptionParser

'''
vcf2bed.py

Simple VCF to BED converter.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    vcf_file = args[0]
    bed_file = args[1]
    bed_out = open(bed_file, 'w')

    for line in open(vcf_file):
    	if not line.startswith('#'):
	    	a = line.split('\t')
	    	chrm = a[0]
	    	pos = int(a[1])
	    	name = a[2]

	    	start = pos - 1
	    	end = start + 1
	    	print('%s\t%d\t%d\t%s' % (chrm,start,end,name), file=bed_out)

    bed_out.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
