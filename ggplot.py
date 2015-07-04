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
def plot(r_script, df_dict, args, df_file=None, print_cmd=False, sep=' '):
    # open temp file
    if df_file == None:
        df_fd, df_file = tempfile.mkstemp()
    else:
        df_fd = None
    df_out = open(df_file, 'w')

    # get headers
    headers = sorted(df_dict.keys())
    print >> df_out, sep.join([str(head) for head in headers])

    # check list lengths
    length = len(df_dict[headers[0]])
    for i in range(1,len(headers)):
        if length != len(df_dict[headers[i]]):
            print >> sys.stderr, 'Lists in dict vary in length.'
            exit(1)

    # print data frame
    for i in range(length):
        print >> df_out, sep.join([str(df_dict[head][i]) for head in headers])
    df_out.close()

    # convert args to one string
    args_str = sep.join([str(a) for a in args])

    # plot in R
    cmd = 'R --slave --args %s %s < %s' % (df_file, args_str, r_script)

    if print_cmd:
        print >> sys.stderr, cmd
        
    subprocess.call(cmd, shell=True)

    # clean
    if df_fd != None:
        os.close(df_fd)
        os.remove(df_file)


################################################################################
# print_df
#
# Just print the given data frame dictionary to the output file given.
################################################################################
def print_df(df_dict, out_file=None):
    # open
    if out_file == None:
        df_fd, df_file = tempfile.mkstemp()
    else:
        df_file = out_file
    df_out = open(df_file, 'w')

    # get headers
    headers = sorted(df_dict.keys())
    print >> df_out, ' '.join([str(head) for head in headers])

    # check list lengths
    length = len(df_dict[headers[0]])
    for i in range(1,len(headers)):
        if length != len(df_dict[headers[i]]):
            print >> sys.stderr, 'Lists in dict vary in length:'
            for j in range(len(headers)):
                print >> sys.stderr, headers[j], len(df_dict[headers[j]])
            exit(1)

    # print data frame
    for i in range(length):
        print >> df_out, ' '.join([str(df_dict[head][i]) for head in headers])
    df_out.close()

    if out_file == None:
        return df_fd, df_file
    else:
        return None
