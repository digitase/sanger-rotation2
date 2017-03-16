'''Function definitions for analyses in ./1_process_clonotype_data.py
'''

import statsmodels.api as sm
import scipy
import numpy as np
import functools
import collections
import math
import re

#
# General
#

def get_clonotype_freq(df, clonotype_field="clonotype", groups=None, prop=False, fillna=True):
    '''Get frequency (or proportion) of the clonotype_field, grouped by groups
    '''
    if groups is None:
        freqs = df.groupby([clonotype_field]).apply(len)
    else:
        freqs = df.groupby([clonotype_field] + groups).apply(len).unstack()
    if fillna:
        freqs = freqs.fillna(0)
    if prop:
        if len(freqs.shape) > 1:
            return freqs.apply(lambda x: x/sum(x)) 
        else:
            return freqs/sum(freqs) 
    else:
        return freqs

def get_naive_rep(df, negate=False):
    '''Filter for the naive repertoire

    negate: If true, return non-naive rep instead.
    '''
    mask = ((df["V-REGION identity %"] == 100)
            & (df["digest"].map(lambda x: "IGHM" in x or "IGHG" in x))
            & (df["cell_type"] == "PBMCs"))
    if negate:
        return df.loc[np.logical_not(mask)]
    else:
        return df[mask]

#
# Specific expansions
#

def arcsin_transform(x):
    ''' Return the arcsine of the square root
    '''
    return(np.arcsin(np.sqrt(x)))

def prop_test_ztest(rep1_freqs, rep2_freqs, alternative="two-sided"):
    '''Z-test for difference between two proportions

    rep1_freqs, rep2_freqs: frequency tables for the two repertoires
    '''
    assert(len(rep1_freqs) == len(rep2_freqs))
    zs = []
    ps = []
    for i in range(len(rep1_freqs)):
        if rep1_freqs[i] or rep2_freqs[i]:
            z, p = sm.stats.proportions_ztest(
                    [rep1_freqs[i], rep2_freqs[i]],
                    [sum(rep1_freqs), sum(rep2_freqs)],
                    alternative=alternative
                )
            zs.append(z)
            ps.append(p)
        # Situation where both counts are 0
        # We append a non significant result instead of an NA,
        # as all clonotypes tested have counts in some sample, 
        # so will be included in the analysis at some point.
        else:
            zs.append(0.0)
            ps.append(1.0)
    return zs, ps

def prop_test_fisher(rep1_freqs, rep2_freqs, fdr_alpha=None):
    '''Fisher's exact test for association between two variables

    rep1_freqs, rep2_freqs: frequency tables for the two repertoires
    fdr_alpha: alpha level to use for FDR correction

    Note: 
        Fisher's exact test assumes that the row and column totals are fixed, or
        "conditioned" [...] When one or both of the row or column totals are
        unconditioned, the Fisher's exact test is not, strictly speaking, exact.
        Instead, it is somewhat conservative [...]
        Source:
            http://www.biostathandbook.com/fishers.html
    '''
    assert(len(rep1_freqs) == len(rep2_freqs))
    ps = []
    for i in range(len(rep1_freqs)):
        # Conduct test by constructing 2x2 contingency table
        cont_table = np.array([[rep1_freqs[i], sum(rep1_freqs) - rep1_freqs[i]], 
                               [rep2_freqs[i], sum(rep2_freqs) - rep2_freqs[i]]])
        ps.append(scipy.stats.fisher_exact(cont_table)[1])
    if fdr_alpha:
        ps = sm.stats.multipletests(ps, method="fdr_bh", alpha=fdr_alpha)
    return ps

def get_mean_ranks(x1, x2):
    '''Get the mean ranks of samples x1 and x2 in the combined sample

    Smallest value = rank 1

    Returns: mean ranks and medians of samples
    '''
    combined = np.concatenate([x1, x2], axis=0)
    ranks = scipy.stats.rankdata(combined)
    ranks1, ranks2 = ranks[:len(x1)], ranks[len(x1):]
    return np.mean(ranks1), np.mean(ranks2), np.median(x1), np.median(x2)

def get_isotype_distribution(clonotype, rep):
    '''Get frequencies of isotypes for a clonotype in a repertoire

    Fractional counts are assigned for clones with multiple isotypes.
    '''
    def list_to_counter(x):
        '''Convert a list to a scaled counter
        '''
        c = collections.Counter(x)
        for i in c:
            c[i] /= len(c)
        return c

    clones = rep.loc[(rep['clonotype'] == clonotype), 'isotypes']
    # Add counters to get clonotype freqs
    return functools.reduce(lambda x, y: x+y, list(clones.apply(list_to_counter)), collections.Counter())

def get_prop_ci(p, n, z=1.96, interval=False):
    '''Get CI for a population proportion'''
    rad = z*math.sqrt(p*(1-p)/n)
    if interval:
        return (p-rad, p+rad)
    else:
        return rad

def natural_key(string_):
    '''Sorting by numeric values in strings

    See http://www.codinghorror.com/blog/archives/001018.html
    ''' 
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

