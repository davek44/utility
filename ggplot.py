#!/usr/bin/env python
from optparse import OptionParser
import os, subprocess, sys, tempfile

################################################################################
# ggplot.py
#
# Make a plot given an R script, dict data frame, and arguments.
################################################################################


################################################################################
# plot
################################################################################
def plot(r_script, df_dict, args, debug=False):
    # open temp file
    if debug:
        df_file = 'data_frame.txt'
    else:
        df_fd, df_file = tempfile.mkstemp()
    df_out = open(df_file, 'w')

    # get headers
    headers = sorted(df_dict.keys())
    print >> df_out, ' '.join([str(head) for head in headers])

    # check list lengths
    length = len(df_dict[headers[0]])
    for i in range(1,len(headers)):
        if length != len(df_dict[headers[i]]):
            print >> sys.stderr, 'Lists in dict vary in length.'
            exit(1)

    # print data frame
    for i in range(length):
        print >> df_out, ' '.join([str(df_dict[head][i]) for head in headers])
    df_out.close()

    # convert args to one string
    args_str = ' '.join([str(a) for a in args])

    # plot in R
    subprocess.call('R --slave --args %s %s < %s' % (df_file, args_str, r_script), shell=True)

    # clean
    if not debug:
        os.close(df_fd)
        os.remove(df_file)
