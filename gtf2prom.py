#!/usr/bin/env python
from optparse import OptionParser
import gff
import os, subprocess

################################################################################
# gtf2prom.py
#
# Produce a GFF file and a fasta file corresponding to the promoter of the
# genes in a gtf fle.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    parser.add_option('-d', dest='downstream', type='int', default=0, help='Downstream promoter length [Default: %default]')
    parser.add_option('-u', dest='upstream', type='int', default=2000, help='Upstream promoter length [Default: %default]')
    parser.add_option('-o', dest='output_pre', default='promoter', help='Output file prefix [Default: %default]')
    (options,args) = parser.parse_args()
    
    if len(args) != 1:
        parser.error('Must provide gtf file')
    else:
        gtf_file = args[0]

    gff.promoters(gtf_file, options.upstream, options.downstream, '%s.gff'%options.output_pre)
    p = subprocess.Popen('gff2fa.py %s.gff > %s.fa' % (options.output_pre,options.output_pre), shell=True)
    os.waitpid(p.pid,0)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
