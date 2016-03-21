#!/usr/bin/env python
from optparse import OptionParser

import matplotlib.pyplot as plt

################################################################################
# plot.py
#
# matplotlib helper methods
################################################################################

#####################################################################
# limits
#
# Determine a nice buffered axis range from a list/array of numbers
#####################################################################
def limits(nums, buf_pct=0.05):
    nmin = min(nums)
    nmax = max(nums)
    spread = nmax-nmin
    buf = buf_pct*spread
    return nmin-buf, nmax+buf

#####################################################################
# scatter
#
# Example scatter plot with some reasonable parameter choices
#####################################################################
def scatter(x, y, pdf, xlabel='', ylabel=''):
    f, ax = plt.subplots()

    # scatter
    plt.scatter(x, y, s=20, alpha=0.8, linewidths=0)

    # x-axis
    xmin, xmax = limits(x)
    plt.xlim(xmin, xmax)
    plt.xlabel(xlabel)
    ax.xaxis.label.set_fontsize(18)
    map(lambda xl: xl.set_fontsize(15), ax.get_xticklabels())

    # y-axis
    ymin, ymax = limits(y)
    plt.ylim(ymin, ymax)
    plt.ylabel(ylabel)
    ax.yaxis.label.set_fontsize(18)
    map(lambda yl: yl.set_fontsize(15), ax.get_yticklabels())

    # save
    plt.savefig(pdf)

    # close
    plt.close()
