#!/usr/bin/env python
import numpy as np
from scipy.stats import rankdata, tiecorrect
from scipy.stats.distributions import norm
import math, os, random, subprocess, tempfile

################################################################################
# stats.py
#
# Common statistical methods.
################################################################################


################################################################################
# entropy
#
# Compute entropy of the given list of numbers.
################################################################################
def entropy(ls):
    return sum([-l*math.log(l) for l in ls if l > 0])


################################################################################
# geo_mean
################################################################################
def geo_mean(ls, log_sum=True, pseudocount=0):
    if log_sum:
        return math.exp(sum([math.log(x+pseudocount) for x in ls])/float(len(ls)))
    else:
        prod = ls[0]
        for x in ls:
            prod *= x
        return math.pow(prod,1.0/len(ls))


################################################################################
# jsd
#
# Jensen-Shannon divergence
################################################################################
def jsd(P, Q):
    if len(P) != len(Q):
        raise ValueError('Distributions P (%d) and Q (%d) are different lengths' % (len(P),len(Q)))

    M = [0.5*P[i]+0.5*Q[i] for i in range(len(P))]

    return 0.5*kld(P,M) + 0.5*kld(Q,M)
    #return entropy(M) - 0.5*entropy(P) - 0.5*entropy(Q)


################################################################################
# kld
#
# Kullback-Leibler divergence
################################################################################
def kld(P, Q):
    if len(P) != len(Q):
        raise ValueError('Distributions P (%d) and Q (%d) are different lengths' % (len(P),len(Q)))
    else:
        invalid = [i for i in range(len(P)) if P[i] > 0 and Q[i] == 0]
        if invalid:
            raise ValueError('Invalid term: P[%d] = %f and Q[%d] = %f' % (invalid[0],P[invalid[0]],invalid[0],Q[invalid[0]]))

    return sum([P[i]*math.log(float(P[i])/Q[i]) for i in range(len(P)) if P[i] > 0])


################################################################################
# mannwhitneyu
#
# My version that returns the z value like ranksums.
################################################################################
def mannwhitneyu(x, y, use_continuity=True):
    """
    Computes the Mann-Whitney rank test on samples x and y.

    Parameters
    ----------
    x, y : array_like
        Array of samples, should be one-dimensional.
    use_continuity : bool, optional
            Whether a continuity correction (1/2.) should be taken into
            account. Default is True.

    Returns
    -------
    u : float
        The Mann-Whitney statistics.
    prob : float
        One-sided p-value assuming a asymptotic normal distribution.

    Notes
    -----
    Use only when the number of observation in each sample is > 20 and
    you have 2 independent samples of ranks. Mann-Whitney U is
    significant if the u-obtained is LESS THAN or equal to the critical
    value of U.

    This test corrects for ties and by default uses a continuity correction.
    The reported p-value is for a one-sided hypothesis, to get the two-sided
    p-value multiply the returned p-value by 2.

    """
    x = np.asarray(x)
    y = np.asarray(y)
    n1 = len(x)
    n2 = len(y)
    ranked = rankdata(np.concatenate((x,y)))
    rankx = ranked[0:n1]       # get the x-ranks
    #ranky = ranked[n1:]        # the rest are y-ranks
    u1 = n1*n2 + (n1*(n1+1))/2.0 - np.sum(rankx,axis=0)  # calc U for x
    u2 = n1*n2 - u1                            # remainder is U for y
    bigu = max(u1,u2)
    smallu = min(u1,u2)
    #T = np.sqrt(tiecorrect(ranked))  # correction factor for tied scores
    T = tiecorrect(ranked)
    if T == 0:
        raise ValueError('All numbers are identical in amannwhitneyu')
    sd = np.sqrt(T*n1*n2*(n1+n2+1)/12.0)

    if use_continuity:
        # normal approximation for prob calc with continuity correction
        z = (bigu-0.5-n1*n2/2.0) / sd
    else:
        z = (bigu-n1*n2/2.0) / sd  # normal approximation for prob calc
    z *= int(u1<u2)-int(u1>u2)
    return z, norm.sf(abs(z))  #(1.0 - zprob(z))


############################################################
# max_i
#
# Find max and return index and value
############################################################
def max_i(lis):
    max_index = 0
    max_val = lis[0]
    for i in range(1,len(lis)):
        if lis[i] > max_val:
            max_index = i
            max_val = lis[i]

    return (max_val,max_index)


############################################################
# min_i
#
# Find min and return index and value
############################################################
def min_i(lis):
    min_index = 0
    min_val = lis[0]
    for i in range(1,len(lis)):
        if lis[i] < min_val:
            min_index = i
            min_val = lis[i]

    return (min_val,min_index)


############################################################
# median
#
# Return the median of a list
############################################################
def median(ls):
    sls = sorted(ls)
    if len(sls) % 2 == 1:
        return sls[(len(sls)+1)/2-1]
    else:
        lower = sls[len(sls)/2-1]
        upper = sls[len(sls)/2]
        return float(lower+upper)/2.0


############################################################
# mean
#
# Return mean of a list
############################################################
def mean(ls):
    if len(ls) == 0:
        return 0
    else:
        return float(sum(ls)) / float(len(ls))

    
############################################################
# mean_sd
#
# Return the mean and sd of a list
############################################################
def mean_sd(ls):
    u = mean(ls)
    dev_sum = 0.0
    for x in ls:
        dev_sum += (x-u)*(x-u)
    return u, math.sqrt(dev_sum / float(len(ls)))


################################################################################
# mi_parmigene
#
# Compute mutual information on continuous arrays using parmigene in R.
################################################################################
def mi_parmigene(array1, array2, debug=False):
    df_dict = {'A':array1, 'B':array2}

    # open temp file
    if debug:
        df_file = 'data_frame.txt'
    else:
        df_fd, df_file = tempfile.mkstemp()
    df_out = open(df_file, 'w')

    # print headers
    print >> df_out, 'A B'
    
    # check list lengths
    length = len(df_dict['A'])
    if length != len(df_dict['B']):
        print >> sys.stderr, 'Lists in dict vary in length.'
        exit(1)

    # print data frame
    for i in range(length):
        print >> df_out, '%s %s' % (str(df_dict['A'][i]), str(df_dict['B'][i]))
    df_out.close()

    # compute in R
    mi = float(subprocess.check_output('R --slave --args %s < %s/mi_parmigene.r' % (df_file,os.environ['RDIR']), shell=True))

    # clean
    if not debug:
        os.close(df_fd)
        os.remove(df_file)

    return mi


################################################################################
# mutual_information
#
# Input given as a discrete probability distribution matrix.
################################################################################
def mutual_information(m):
    n = m.shape[0]

    # compute single variable distributions
    px = [0]*n
    py = [0]*n
    for i in range(n):
        px[i] = sum(m[i,:])
        py[i] = sum(m[:,i])

    # sum mutual information
    mi = 0.0
    for i in range(n):
        for j in range(n):
            if m[i,j]:
                mi += m[i,j]*math.log(float(m[i,j])/(px[i]*py[j]))

    return mi


################################################################################
# normalize
#
# To sum to 1.
################################################################################
def normalize(ls):
    ls_sum = float(sum(ls))
    return [l/ls_sum for l in ls]


############################################################
# quantile
#
# Return the value at the quantile given.
############################################################
def quantile(ls, q):
    sls = sorted(ls)

    if type(q) == list:
        qval = []
        for j in range(len(q)):
            qi = int(len(sls)*q[j] + 0.5)
            qval.append(sls[qi])
    else:
        qi = int(len(sls)*q + 0.5)
        qval = sls[qi]

    return qval


############################################################
# sd
#
# Return the standard deviation of a list
############################################################
def sd(ls):
    return math.sqrt(variance(ls))


################################################################################
# sample_probs
#
# Sample from a list of items according to given probabilities
################################################################################
def sample_probs(items, probs, count=1):
    # compute cumulative probabilities
    cum_probs = [0]*len(probs)
    cum_probs[0] = probs[0]
    for p in range(1,len(probs)):
        cum_probs[p] = cum_probs[p-1] + probs[p]

    # sample
    samples = ['']*count
    for c in range(count):
        # get a random number
        r = random.random()
        p = 0
        while p < len(cum_probs) and r > cum_probs[p]:
            p += 1
        samples[c] = items[p]

    return samples


############################################################
# variance
#
# Return the variance of a list.
############################################################
def variance(ls):
    u = mean(ls)
    dev_sum = 0.0
    for x in ls:
        dev_sum += (x-u)*(x-u)
    return dev_sum / float(len(ls)-1)
