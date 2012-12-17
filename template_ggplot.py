#!/usr/bin/env python
from optparse import OptionParser
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.lib.ggplot2 as ggplot2

grdevices = importr('grDevices')

################################################################################
# name
#
#
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    #parser.add_option()
    (options,args) = parser.parse_args()

    # construct data frame
    df = ro.DataFrame({'a':a, 'b':b})

    # construct full plot
    gp = ggplot2.ggplot(df) + \
        ggplot2.aes_string(x='a', y='b') + \
        ggplot2.geom_point() + \
        ggplot2.scale_x_continuous('a') + \
        ggplot2.scale_y_continuous('b')

    # plot to file
    grdevices.pdf(file='temp.pdf')
    gp.plot()
    grdevices.dev_off()

    

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
