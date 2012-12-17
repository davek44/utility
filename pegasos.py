#!/usr/bin/env python
from optparse import OptionParser
from scipy import zeros
from scipy.stats import binom
from numpy import dot
import math, os, pdb, random, shutil, sys
import dna, stats, util

################################################################################
# pegasos.py
#
# Code to aid use of the Pegasos SVM package.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <input_file>'
    parser = OptionParser(usage)
    parser.add_option('-k', dest='k_fold', type='int', default=10, help='Number of folds to use for cross-validation [Default: %default]')
    parser.add_option('--lambda_min', dest='lambda_min', type='float', default=.01, help='Minimum -lambda value to attempt [Default: %default]')
    parser.add_option('--lambda_max', dest='lambda_max', type='float', default=10.0, help='Maximum -lambda value to attempt [Default: %default]')
    parser.add_option('--lambda_mult', dest='lambda_mult', type='float', default=2.0, help='Multiplier for next -lambda value to attempt [Default: %default]')
    parser.add_option('-l', dest='lesser_kmers', action='store_true', default=False, help='Use all kmers of length less than and equal to that given by -k [Default: %default]')
    #parser.add_option('-m', dest='model_file', help='File to output model to')
    parser.add_option('-p', dest='parallel', type='int', default=4, help='Number of parallel threads to run [Default: %default]')
    parser.add_option('-r', dest='replicates', type='int', default=1, help='Number of times to repeat the optimization for each fold [Default: %default]')
    parser.add_option('-w', dest='weights', action='store_true', default=False, help='Print a summary of the weight vectors')
    (options,args) = parser.parse_args()

    if len(args) != 1:
        parser.error('Must provide input file')
    else:
        input_file = args[0]
    input_base = os.path.splitext(input_file)[0]

    if options.weights:
        summarize_weights(input_base, options)
        exit()

    # determine % of positive examples
    input_pos, input_total = positive_percent(input_file)
    f1_base = input_pos / float(input_total) # trust me, it works

    for r in range(options.replicates):
        rep_dir = '%s_rep%d' % (input_base,r)
        if os.path.isdir(rep_dir):
            shutil.rmtree(rep_dir)
        os.mkdir(rep_dir)
        os.chdir(rep_dir)

        # divide data into folds
        divide_data('../'+input_file, options.k_fold)

        # collect pegasos commands
        cmds = []
        peg_lambda = options.lambda_min
        while peg_lambda <= options.lambda_max:
            # run on each fold
            for f in range(options.k_fold):
                cmds.append('pegasos -lambda %f -modelFile fold%d/train_%.1e.mod fold%d/train.dat &> /dev/null' % (peg_lambda, f, peg_lambda, f))

            # increase lambda
            peg_lambda *= options.lambda_mult

        # exceute pegasos commands
        util.exec_par(cmds, options.parallel)

        # start to clean up space
        for f in range(options.k_fold):
            os.remove('fold%d/train.dat' % f)

        os.chdir('..')


    # collect results
    peg_lambda = options.lambda_min
    while peg_lambda <= options.lambda_max:
        recalls = []
        precisions = []
        failed = False

        for r in range(options.replicates):
            if not failed:
                outcomes = {'tp':0, 'fp':0, 'fn':0}

                # collect each fold
                for f in range(options.k_fold):
                    if not compute_accuracy(outcomes, '%s_rep%d/fold%d'%(input_base,r,f), peg_lambda):
                        failed = True
                        break

                # save
                if not failed:
                    recalls.append(float(outcomes['tp'])/(outcomes['tp']+outcomes['fn']))
                    precisions.append(float(outcomes['tp'])/(outcomes['tp']+outcomes['fp']))

        # summarize and print
        if failed:
            print '%.1e %8s %7s %8s %7s %8s %8s' % (peg_lambda, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA')
        else:
            recall, rsd = stats.mean_sd(recalls)
            rsd /= math.sqrt(len(recalls))
            precision, psd = stats.mean_sd(precisions)
            psd /= math.sqrt(len(precisions))

            #null_p = 1.0-binom.cdf(int(recall*input_total+0.5)-1, int(recall*input_total/precision + 0.5), float(input_pos)/input_total)

            f1 = 2*recall*precision/(recall+precision)            

            #print '%.1e %8.3f %6.3f %8.3f %6.3f %8.3f %8.3f %8.1e' % (peg_lambda, recall, rsd, precision, psd, f1, (f1-f1_base), null_p)
            print '%.1e %8.4f %7.4f %8.4f %7.4f %8.4f %8.4f' % (peg_lambda, recall, rsd, precision, psd, f1, (f1-f1_base))

        peg_lambda *= options.lambda_mult



################################################################################
# compute_accuracy
#
# Process the test file defined by fold using the model file defined by fold
# peg_lambda and compute accuracy statistics.
################################################################################
def compute_accuracy(outcomes, fold_dir, peg_lambda, b=0):
    # read model
    model_line = open('%s/train_%.1e.mod' % (fold_dir,peg_lambda)).readline()
    w_entries = model_line.split()
    w = zeros(len(w_entries))
    for w_entry in w_entries:
        i, wi = w_entry.split(':')
        w[int(i)-1] = float(wi)

    if len([wi for wi in w if math.isinf(wi) or math.isnan(wi)]):
        # verify model has no infinities
        return False
    
    else:
        # process test set
        tp = 0
        fp = 0
        fn = 0
        for line in open('%s/test.dat' % fold_dir):
            a = line.split()

            # read test entry
            clas = int(a[0])
            x = zeros(len(w))
            for x_entry in a[1:]:
                i, xi = x_entry.split(':')
                x[int(i)-1] = float(xi)

            # predict class and compare
            wx = dot(w,x)
            if wx > b:
                if clas == 1:
                    tp += 1
                else:
                    fp += 1
            elif clas == 1:
                fn += 1

        # update counts
        outcomes['tp'] += tp
        outcomes['fp'] += fp
        outcomes['fn'] += fn

        return True


################################################################################
# divide_data
#
# Divide the input data into folds for training and testing.
################################################################################
def divide_data(input_file, k_fold):
    # read data
    input_lines = open(input_file).readlines()

    test_num = int(len(input_lines)*1.0/k_fold)

    for f in range(k_fold):
        # create directory
        if os.path.isdir('fold%d' % f):
            shutil.rmtree('fold%d' % f)
        os.mkdir('fold%d' % f)

        # sample indexes for testing
        test_indexes = sorted(random.sample(range(len(input_lines)), test_num))

        # divide input between test and training
        ti = 0
        test_out = open('fold%d/test.dat' % f, 'w')
        train_out = open('fold%d/train.dat' % f, 'w')
        for i in range(len(input_lines)):
            if ti < len(test_indexes) and i == test_indexes[ti]:
                print >> test_out, input_lines[i],
                ti += 1
            else:
                print >> train_out, input_lines[i],

        test_out.close()
        train_out.close()


################################################################################
# positive_percent
#
# Return the proportion of positive examples in the data.
################################################################################
def positive_percent(input_file):
    pos = 0
    total = 0
    for line in open(input_file):
        a = line.split()
        if int(a[0]) > 0:
            pos += 1
        total += 1
    return pos, total
   

################################################################################
# summarize_weights
#
# Print a summary table about the learned weights.
################################################################################
def summarize_weights(input_base, options):
    peg_lambda = options.lambda_min
    while peg_lambda <= options.lambda_max:
        # read weights from training models
        weights = []
        for r in range(options.replicates):
            for f in range(options.k_fold):
                line = open('%s_rep%d/fold%d/train_%.1e.mod' % (input_base, r, f, peg_lambda)).readline()
                a = line.split()
                for i_w in a:
                    (i,w) = i_w.split(':')
                    while int(i)-1 >= len(weights):
                        weights.append([])
                    weights[int(i)-1].append(float(w))

        # determine kmers
        if options.lesser_kmers:
            min_k = 1
            max_k = int(math.log(1+3*(len(weights)+1), 4) - 1 + 0.5) # 0.5 for rounding
        else:
            min_k = int(math.log(len(weights), 4) + 0.5)
            max_k = min_k

        kmers = []
        for k in range(min_k,max_k+1):
            for i in range(int(math.pow(4,k))):
                kmers.append(dna.int2kmer(k,i))
        kmers.sort()

        # summarize and print
        for i in range(len(weights)):
            if options.lambda_min == options.lambda_max:
                print '%-7s %8.4f %8.4f' % (kmers[i], stats.mean(weights[i]), stats.sd(weights[i]))
            else:
                print '%-7.2e %-7s %8.4f %8.4f' % (peg_lambda, kmers[i], stats.mean(weights[i]), stats.sd(weights[i])/math.sqrt(len(weights[i])))

        peg_lambda *= options.lambda_mult
     
################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
