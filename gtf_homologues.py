#!/usr/bin/env python
from optparse import OptionParser
import subprocess
import gff

################################################################################
# gtf_homologues.py
#
# Make a table describing candidate homologue genes as determined by a
# transmap from one genome to another.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <chain_file> <net_file> <gtf_from> <gtf_to>'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 4:
        parser.error('Must provide chain file and two GTF files')
    else:
        chain_file = args[0]
        net_file = args[1]
        gtf_from = args[2]
        gtf_to = args[3]

    # transmap to new genome
    subprocess.call('chain_map.py -k gene_id -n %s %s %s > from_map.gtf' % (net_file,chain_file,gtf_from), shell=True)

    # intersect w/ gtf_to
    homologues = {}
    p = subprocess.Popen('intersectBed -wo -s -a from_map.gtf -b %s' % gtf_to, shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        a = line.split('\t')
        
        kv_to = gff.gtf_kv(a[17])

        gid_from = a[8].split(';')[1].strip()
        gid_to = kv_to['gene_id']

        homologues.setdefault(gid_from,set()).add(gid_to)
    p.communicate()

    # find all genes
    genes = set()
    for line in open(gtf_from):
        a = line.split('\t')
        genes.add(gff.gtf_kv(a[8])['gene_id'])

    # print table
    for g in genes:
        print '%s\t%s' % (g,' '.join(homologues.get(g,[])))


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
