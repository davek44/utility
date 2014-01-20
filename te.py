#!/usr/bin/env python
from optparse import OptionParser
import glob, os, sys
import gff

################################################################################
# te.py
#
# Methods to work with transposable element annotations.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()
    

################################################################################
# map_dfam_repeat
#
# Return a dict mapping DFAM repeats to RepeatMasker repeats.
################################################################################
def map_dfam_repeat():
    repeats = set()
    for line in open('%s/hg19.fa.out.tp.gff' % os.environ['MASK']):
        a = line.split('\t')
        kv = gff.gtf_kv(a[8])
        repeats.add(kv['repeat'])

    dfam_repeat = {}
    for repeat in repeats:
        dfam_tes = map_rm_dfam(repeat, quiet=True)
        for dfam_te in dfam_tes:
            dfam_repeat[dfam_te] = repeat

    return dfam_repeat


################################################################################
# map_rm_dfam
#
# Map a RepeatMasker name to a DFAM name.
################################################################################
def map_rm_dfam(repeat, quiet=False):
    if os.path.isfile('%s/hmms/%s.hmm' % (os.environ['DFAM'],repeat)):
        dfam_reps = [repeat]
    elif os.path.isfile('%s/hmms/%sv.hmm' % (os.environ['DFAM'],repeat)):
        dfam_reps = [repeat+'v']
    else:
        hmm_files = glob.glob('%s/hmms/%s_*.hmm' % (os.environ['DFAM'],repeat))

        # if no hits
        if len(hmm_files) == 0:
            # try removing "-int"
            if repeat[-4:] == '-int' and os.path.isfile('%s/hmms/%s.hmm' % (os.environ['DFAM'],repeat[:-4])):
                dfam_reps = [repeat[:-4]]
            else:
                # missing
                if not quiet:
                    print >> sys.stderr, 'Missing DFAM name for %s' % repeat
                dfam_reps = []

        # if hits
        else:
            # grab em
            dfam_reps = []
            for i in range(len(hmm_files)):
                start = hmm_files[i].rfind('/')+1
                end = hmm_files[i].rfind('.hmm')
                dfam_reps.append(hmm_files[i][start:end])

    return dfam_reps


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
