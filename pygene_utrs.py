#!/usr/bin/env python
from optparse import OptionParser

import pygene

'''
Name

Description...
'''

################################################################################
# main
################################################################################
def main():
  usage = 'usage: %prog [options] <in_gtf_file> <out_gtf_file>'
  parser = OptionParser(usage)
  (options,args) = parser.parse_args()

  if len(args) != 2:
  	parser.error('Must provide input and output GTFs')
  else:
  	in_gtf_file = args[0]
  	out_gtf_file = args[1]

  gtf = pygene.GTF(in_gtf_file)

  with open(out_gtf_file, 'w') as out_gtf_open:
    gtf.write_gtf(out_gtf_file, write_cds=True, write_utrs=True)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
  main()
