#!/usr/bin/env python
from optparse import OptionParser
import numpy as np

'''
sciseq_collision.py

Estimate the collision rate for a set of sci-seq barcode parameters.
'''

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] arg'
    parser = OptionParser(usage)
    parser.add_option('-b', dest='barcode1',
            default=None, type='int',
            help='Number of barcodes introduced in the first RT stage')
    parser.add_option('-c', dest='cells',
            default=None, type='int',
            help='Number of cells sorted per well in the second PCR stage')
    parser.add_option('-n', dest='num_samples',
            default=10000, type='int',
            help='Number of simulation samples [Default: %default]')
    (options,args) = parser.parse_args()

    collisions = 0

    for i in range(options.num_samples):
        cell_barcodes = np.random.randint(0, options.barcode1, size=options.cells)
        unique_cell_barcodes = len(set(cell_barcodes))
        cell_collisions = options.cells - unique_cell_barcodes
        collisions += cell_collisions

    collision_rate = collisions / (options.num_samples*options.cells)
    print('Collision rate: %.4f' % collision_rate)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
