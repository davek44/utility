#!/usr/bin/env python
from optparse import OptionParser

import numpy as np

################################################################################
# meme2possum.py
#
# Convert a file of MEME PWMs to Possum's input format.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <meme_file> <possum_file>'
    parser = OptionParser(usage)
    # parser.add_option()
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide input MEME file and output Possum file')
    else:
        meme_file = args[0]
        possum_file = args[1]

    #######################################################
    # input MEME motifs
    #######################################################
    motif_pwms = {}
    in_motif = False
    for line in open(meme_file):
        if line.startswith('MOTIF'):
            motif_id = line.split()[1]
            in_motif = True
            pwm_cols = []
        elif in_motif:
            if line.startswith('letter-probability matrix'):
                pass
            elif line.strip() == '':
                in_motif = False
                motif_pwms[motif_id] = np.array(pwm_cols)
            else:
                pwm_cols.append([float(p) for p in line.split()])

    if in_motif:
        motif_pwms[motif_id] = np.array(pwm_cols)

    #######################################################
    # output Possum
    #######################################################
    possum_out = open(possum_file, 'w')
    print >> possum_out, 'BEGIN GROUP'

    for motif_id in motif_pwms:
        mpwm = motif_pwms[motif_id]
        motif_len = mpwm.shape[0]

        print >> possum_out, 'BEGIN FLOAT'
        print >> possum_out, 'ID %s' % motif_id
        print >> possum_out, 'AP DNA'
        print >> possum_out, 'LE %d' % motif_len
        for ci in range(motif_len):
            print >> possum_out, 'MA %s' % ' '.join([str(n) for n in mpwm[ci]])
        print >> possum_out, 'END'
        print >> possum_out, ''

    print >> possum_out, 'END'

    possum_out.close()


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
