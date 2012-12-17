#!/usr/bin/env python
from optparse import OptionParser
import gff
import os, pdb, subprocess

################################################################################
# merge_gff_spliced.py
#
# Similar to mergeBed but additional functionality to handle spliced
# features.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <gff file>'
    parser = OptionParser(usage)
    parser.add_option('-k', dest='key', help='Key used to link feature entries [Default: %default]')
    parser.add_option('-p', dest='pct_t', default=.9, type='float', help='Percentage coverage threshold to merge [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide gff file')
    else:
        gff_file = args[0]
        gff_pre = os.path.splitext(gff_file)[0]

    # get features
    feature_lengths = {}
    for line in open(gff_file):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        key = a[8]
        if options.key:
            key = gff.gtf_kv(a[8])[options.key]

        feature_lengths[key] = feature_lengths.get(key,0) + int(a[4])-int(a[3])+1
        
    # intersect features
    p = subprocess.Popen('intersectBed -s -wo -a %s -b %s > %s_self.gff' % (gff_file,gff_file,gff_pre), shell=True)
    os.waitpid(p.pid,0)

    # hash overlap bp
    overlap_bp = {}
    for line in open('%s_self.gff' % gff_pre):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        key1 = a[8]
        if options.key:
            key1 = gff.gtf_kv(a[8])[options.key]
        key2 = a[17]
        if options.key:
            key2 = gff.gtf_kv(a[17])[options.key]

        if key1 < key2: # just in one direction
            if key1 not in overlap_bp:
                overlap_bp[key1] = {}
            overlap_bp[key1][key2] = overlap_bp[key1].get(key2,0) + int(a[-1])

    # create list of % overlaps
    overlap_pcts = []
    for key1 in overlap_bp:
        for key2 in overlap_bp[key1]:
            pct1 = float(overlap_bp[key1][key2])/feature_lengths[key1]
            pct2 = float(overlap_bp[key1][key2])/feature_lengths[key2]
            max_pct = max(pct1,pct2)
            if max_pct > options.pct_t:
                overlap_pcts.append((max_pct,key1,key2))

    # greedily merge
    overlap_pcts.sort(reverse=True)
    for i in range(len(overlap_pcts)):
        (pct,key1,key2) = overlap_pcts[i]

        if key1 != key2: # means list had a-b, a-c, b-c
            # decided which to delete
            if feature_lengths[key1] > feature_lengths[key2]:
                del_key = key2
                replace_key = key1
            else:
                del_key = key1
                replace_key = key2

            # delete
            del feature_lengths[del_key]

            # replace remainder
            for j in range(i+1,len(overlap_pcts)):
                if del_key == overlap_pcts[j][1]:
                    overlap_pcts[j] = (overlap_pcts[j][0],replace_key,overlap_pcts[j][2])
                if del_key == overlap_pcts[j][2]:
                    overlap_pcts[j] = (overlap_pcts[j][0],overlap_pcts[j][1],replace_key)

    # print un-deleted
    for line in open(gff_file):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        key = a[8]
        if options.key:
            key = gff.gtf_kv(a[8])[options.key]

        if key in feature_lengths:
            print line,


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
