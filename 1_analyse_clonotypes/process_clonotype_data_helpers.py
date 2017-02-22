'''Function definitions for analyses in ./1_process_clonotype_data.py
'''

import statsmodels.api as sm
import scipy
import numpy as np

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

#
# Specific expansions
#

def arcsin_transform(x):
    ''' Return the arcsine of the square root
    '''
    return(np.arcsin(np.sqrt(x)))

def prop_test_ztest(rep1_freqs, rep2_freqs, alternative="two-sided", fdr_alpha=None):
    '''Z-test for difference between two proportions

    rep1_freqs, rep2_freqs: frequency tables for the two repertoires
    fdr_alpha: alpha level to use for FDR correction
    '''
    assert(len(rep1_freqs) == len(rep2_freqs))
    ps = []
    for i in range(len(rep1_freqs)):
        if rep1_freqs[i] or rep2_freqs[i]:
            p = sm.stats.proportions_ztest(
                    [rep1_freqs[i], rep2_freqs[i]],
                    [sum(rep1_freqs), sum(rep2_freqs)],
                    alternative=alternative
                )[1]
            ps.append(p)
        # Situation where both counts are 0
        else:
            ps.append(1.0)
    if fdr_alpha:
        ps = sm.stats.multipletests(ps, method="fdr_bh", alpha=fdr_alpha)
    return ps

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

