#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess, tempfile

################################################################################
# fire_map.py
#
# Map motifs discovered by FIRE and specified in the .signif file against the
# given possum index.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <possum_index> <fire_signif>'
    parser = OptionParser(usage)
    parser.add_option('-r', dest='robust', type='int', default=8, help='Minimum motif robustness score (out of 10) [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide possum index and FIRE .signif')
    else:
        possum_index = args[0]
        fire_signif = args[1]

    # convert the regular expressions to uniform PWMs
    motifs_pwm = read_motif_pwms(fire_signif, options.robust)

    # run possum
    possum_fd, possum_file = run_possum(motifs_pwm, possum_index)

    # convert output to gff
    possum2gff(possum_file)

    os.close(possum_fd)
    os.remove(possum_file)


################################################################################
# read_motif_pwms
#
# Read in the motifs from FIRE output as regular expression motifs and convert
# them to PWMs.
################################################################################
def read_motif_pwms(fire_signif, min_robust):
    motifs_re = []
    for line in open(fire_signif):
        a = line.split()
        robust = int(a[4])
        if robust >= min_robust:
            motifs_re.append(a[0])

    motifs_pwm = {}
    for mre in motifs_re:
        mpwm = []

        i = 0
        while i < len(mre):
            if mre[i] != '.':
                mpwm.append({'A':0, 'C':0, 'G':0, 'T':0})

                if mre[i] == '[':
                    i += 1
                    while mre[i] != ']':
                        mpwm[-1][mre[i]] += 1
                        i += 1
                else:
                    mpwm[-1][mre[i]] += 1
            i += 1

        motifs_pwm[mre] = mpwm
    
    return motifs_pwm


################################################################################
# run_possum
################################################################################
def run_possum(motifs_pwm, possum_index):
    ############################################
    # print pwm's for possum
    pwm_fd, pwm_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
    pwm_out = open(pwm_file, 'w')

    print >> pwm_out, 'BEGIN GROUP'
    for motif_re in motifs_pwm:
        motif_pwm = motifs_pwm[motif_re]

        print >> pwm_out, 'BEGIN INT'
        print >> pwm_out, 'ID %s' % motif_re
        print >> pwm_out, 'AP DNA'
        print >> pwm_out, 'LE %d' % len(motif_pwm)
        for i in range(len(motif_pwm)):
            line = 'MA'
            for nt in ['A','C','G','T']:
                line += ' %d' % motif_pwm[i][nt]
            print >> pwm_out, line
        print >> pwm_out, 'END'
    print >> pwm_out, 'END'
    pwm_out.close()

    ############################################
    # run possum
    possum_fd, possum_file = tempfile.mkstemp(dir='%s/research/scratch/temp' % os.environ['HOME'])
    subprocess.call('possumsearch -pr %s -db %s -freq %s_freqs.txt -lazy -esa -pval 1e-3 -fn -rc -format tabs > %s' % (pwm_file,possum_index,possum_index,possum_file), shell=True)

    # clean
    os.close(pwm_fd)
    os.remove(pwm_file)

    return possum_fd, possum_file


################################################################################
# possum2gff
#
# Convert possum output to gff.
################################################################################
def possum2gff(possum_file):
    for line in open(possum_file):
        a = line.split('\t')
        a[-1] = a[-1].rstrip()

        motif_re = a[0]
        try:
            start = int(a[5])-1
        except:
            print 'ERROR: cant extract start'
            print a
            exit(1)
        end = start+int(a[6])-1
        fnrc = a[7]
        try:
            seq_id = a[16][:a[16].find('.')]
        except:
            print 'ERROR: cant extract seq_id'
            print a
            exit(1)

        if fnrc == 'fn':
            strand = '+'
        else:
            strand = '-'

        out_a = [seq_id, 'possum', 'motif', str(start), str(end), '.', strand, '.', motif_re]
        print '\t'.join(out_a)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
