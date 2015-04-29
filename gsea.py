#!/usr/bin/env python
from optparse import OptionParser
import glob, os, subprocess, sys

################################################################################
# gsea.py
#
# Helper script to run GSEA from CuffDiff output.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <diff>'    
    parser = OptionParser(usage)
    parser.add_option('-c', dest='scheme', default='weighted', help='weighted or classic [Default: %default]')
    parser.add_option('-o', dest='out_dir', default='.')
    parser.add_option('-s', dest='gene_set', default='go', help='Gene sets [Default: %default]')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide .diff')
    else:
        diff_file = args[0]

    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    # choose chip
    gsea_jar = glob.glob('%s/gsea*.jar' % os.environ['GSEA'])[0]
    chip_file = '%s/GENE_SYMBOL.chip' % os.environ['GSEA']

    # choose sets
    if options.gene_set.lower() in ['c5', 'go']:
    	sets_file = '%s/sets/c5.all.v5.0.symbols.gmt' % os.environ['GSEA']
    else:
    	print >> sys.stderr, 'Unrecognized gene set: %s' % options.gene_set
    	exit(1)

    # make rank files
    rank_cmd = 'gsea_rnk.py -o %s %s' % (options.out_dir, diff_file)
    subprocess.call(rank_cmd, shell=True)

    for rank_file in glob.glob('%s/*.rnk' % options.out_dir):
        rank_name = rank_file.split('/')[-1][:-4]

        # run GSEA
        gsea_cmd = 'java -cp %s -Xmx4000m xtools.gsea.GseaPreranked -gmx %s -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk %s -scoring_scheme %s -rpt_label %s -chip %s -include_only_symbols true -make_sets true -plot_top_x 50 -rnd_seed timestamp -set_max 1000 -set_min 10 -zip_report false -out %s -gui false' % (gsea_jar, sets_file, rank_file, options.scheme, rank_name, chip_file, options.out_dir)
        subprocess.call(gsea_cmd, shell=True)

        # consider making a new excel file from the gsea_report_for_na_neg*.xls file
        #  where I strip out the redundant col 1 and stupid col 2.


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
