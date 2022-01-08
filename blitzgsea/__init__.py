import random
import time
import numpy as np
import pandas as pd
from loess.loess_1d import loess_1d
from collections import Counter
from scipy import interpolate

from scipy.stats import kstest
from matplotlib import pyplot as plt
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
import warnings
import multiprocessing

from mpmath import mp
from mpmath import mpf
from mpsci.distributions.normal import invcdf
from mpsci.distributions.gamma import cdf as gammacdf

from scipy.stats import norm
from scipy.stats import gamma

import blitzgsea.enrichr
import blitzgsea.plot

from importlib import reload
reload(blitzgsea.enrichr)
reload(blitzgsea.plot)

mp.dps = 1000
mp.prec = 1000

def strip_gene_set(signature, signature_genes, gene_set):
    return [x for x in gene_set if x in signature_genes]

def enrichment_score(signature, abs_signature, signature_map, gene_set):
    gs = set(gene_set)
    hits = [signature_map[x] for x in gs if x in signature_map.keys()]
    hit_indicator = np.zeros(signature.shape[0])
    hit_indicator[hits] = 1
    no_hit_indicator = 1 - hit_indicator
    number_hits = len(hits)
    number_miss = signature.shape[0] - number_hits
    sum_hit_scores = np.sum(abs_signature.iloc[hits])
    norm_hit =  float(1.0/sum_hit_scores)
    norm_no_hit = float(1.0/number_miss)
    running_sum = np.cumsum(hit_indicator * abs_signature.iloc[:,0] * norm_hit - no_hit_indicator * norm_no_hit)
    nn = np.argmax(np.abs(running_sum))
    es = running_sum[nn]
    return running_sum, es

def get_leading_edge(runningsum, signature, gene_set, signature_map):
    gs = set(gene_set)
    hits = [signature_map[x] for x in gs if x in signature_map.keys()]
    rmax = np.argmax(runningsum)
    rmin = np.argmin(runningsum)
    lgenes = []
    if runningsum[rmax] > np.abs(runningsum[rmin]):
        lgenes = set(hits).intersection(set(range(rmax)))
    else:
        lgenes = set(hits).intersection(set(range(rmin, len(runningsum))))
    return ",".join(signature.index[list(lgenes)])

def get_peak_size(signature, abs_signature, signature_map, size, permutations, seed):
    es = []
    random.seed(seed)
    for _ in range(permutations):
        rgenes = random.sample(list(signature.index), size)
        es.append(enrichment_score(signature, abs_signature, signature_map, rgenes)[1])
    return es

def loess_interpolation(x, y, frac=0.5):
    yl = np.array(y)
    xout, yout, wout = loess_1d(x, yl, frac=frac)
    return interpolate.interp1d(xout, yout)

def estimate_parameters(signature, abs_signature, signature_map, library, permutations: int=2000, symmetric: bool=False, calibration_anchors: int=20, plotting: bool=False, processes=4, verbose=False, seed: int=0):
    ll = []
    for key in library.keys():
        ll.append(len(library[key]))
    cc = Counter(ll)
    set_sizes = pd.DataFrame(list(cc.items()),columns = ['set_size','count']).sort_values("set_size")
    set_sizes["cumsum"] = np.cumsum(set_sizes.iloc[:,1])
    
    ll = [len(library[l]) for l in library]
    nn = np.percentile(ll, q=np.linspace(2, 100, calibration_anchors))
    x = sorted(list(set(np.append([1,4,6, np.max(ll)], nn).astype("int"))))

    jobs = processes
    with multiprocessing.Pool(jobs) as pool:
        args = [(signature, abs_signature, signature_map, xx, permutations, symmetric, seed+xx) for xx in x]
        results = list(tqdm(pool.imap(estimate_anchor_star, args), desc="Calibration", total=len(args)))

    alpha_pos = []
    beta_pos = []
    ks_pos = []
    alpha_neg = []
    beta_neg = []
    ks_neg = []
    pos_ratio = []

    for res in results:
        f_alpha_pos, f_beta_pos, f_ks_pos, f_alpha_neg, f_beta_neg, f_ks_neg, f_pos_ratio = res
        alpha_pos.append(f_alpha_pos)
        beta_pos.append(f_beta_pos)
        ks_pos.append(f_ks_pos)
        alpha_neg.append(f_alpha_neg)
        beta_neg.append(f_beta_neg)
        ks_neg.append(f_ks_neg)
        pos_ratio.append(f_pos_ratio)

    if np.max(pos_ratio) > 1.5 and verbose:
        print('Significant unbalance between positive and negative enrichment scores detected. Signature values are not centered close to 0.')

    x = np.array(x, dtype=float)
    
    f_alpha_pos = loess_interpolation(x, alpha_pos)
    f_beta_pos = loess_interpolation(x, beta_pos, frac=0.2)
    
    f_alpha_neg = loess_interpolation(x, alpha_neg)
    f_beta_neg = loess_interpolation(x, beta_neg, frac=0.2)

    # fix issue with numeric instability
    pos_ratio = pos_ratio - np.abs(0.0001*np.random.randn(len(pos_ratio)))
    f_pos_ratio = loess_interpolation(x, pos_ratio)
    
    if plotting:
        xx = np.linspace(min(x), max(x), 1000)
        
        plt.figure(0)
        
        yy = f_alpha_pos(xx)
        plt.figure(1)
        plt.plot(xx, yy, '--', lw=3)
        plt.plot(x, alpha_pos, 'ko')
        
        yy = f_alpha_neg(xx)
        plt.figure(1)
        plt.plot(xx, yy, '--', lw=3, c="orange")
        plt.plot(x, alpha_neg, 'o', c="coral")

        yy = f_beta_pos(xx)
        plt.figure(1)
        plt.plot(xx, yy, '--', lw=3)
        plt.plot(x, beta_pos, 'ko')
        
        yy = f_beta_neg(xx)
        plt.figure(1)
        plt.plot(xx, yy, '--', lw=3, c="orange")
        plt.plot(x, beta_neg, 'o', c="coral")
        
        yy = f_pos_ratio(xx)
        plt.figure(2)
        plt.plot(xx, yy, lw=3)
        plt.plot(x, pos_ratio, 'o', c="black")
        
    return f_alpha_pos, f_beta_pos, f_pos_ratio, np.mean(ks_pos), np.mean(ks_neg)

def estimate_anchor_star(args):
    return estimate_anchor(*args)

def estimate_anchor(signature, abs_signature, signature_map, set_size, permutations, symmetric, seed):
    es = np.array(get_peak_size(signature, abs_signature, signature_map, set_size, permutations, seed))
    pos = [x for x in es if x > 0]
    neg = [x for x in es if x < 0]
    if (len(neg) < 250 or len(pos) < 250) and not symmetric:
        symmetric = True
    if symmetric:
        aes = np.abs(es)
        fit_alpha, fit_loc, fit_beta = gamma.fit(aes, floc=0)
        ks_pos = kstest(aes, 'gamma', args=(fit_alpha, fit_loc, fit_beta))[1]
        ks_neg = kstest(aes, 'gamma', args=(fit_alpha, fit_loc, fit_beta))[1]

        alpha_pos = fit_alpha
        beta_pos = fit_beta
        
        alpha_neg = fit_alpha
        beta_neg = fit_beta
    else:
        fit_alpha, fit_loc, fit_beta = gamma.fit(pos, floc=0)
        ks_pos = kstest(pos, 'gamma', args=(fit_alpha, fit_loc, fit_beta))[1]
        alpha_pos = fit_alpha
        beta_pos = fit_beta
        
        fit_alpha, fit_loc, fit_beta = gamma.fit(-np.array(neg), floc=0)
        ks_neg = kstest(-np.array(neg), 'gamma', args=(fit_alpha, fit_loc, fit_beta))[1]
        alpha_neg = fit_alpha
        beta_neg = fit_beta
    pos_ratio = len(pos)/(len(pos)+len(neg))

    return alpha_pos, beta_pos, ks_pos, alpha_neg, beta_neg, ks_neg, pos_ratio

def probability_star(args):
    return probability(*args)
    
def probability(signature, abs_signature, signature_map, gene_set, f_alpha_pos, f_beta_pos, f_pos_ratio):
    gsize = len(gene_set)
    
    rs, es = enrichment_score(signature, abs_signature, signature_map, gene_set)
    #legenes = get_leading_edge(rs, signature, gene_set, signature_map)

    #pos_alpha = f_alpha_pos(gsize)
    #pos_beta = f_beta_pos(gsize)
    #pos_ratio = f_pos_ratio(gsize)

    nes = 0
    pval = 1
    legenes = ""

    return gsize, es, nes, pval, legenes

def gsea(signature, library, permutations: int=2000, anchors: int=20, min_size: int=5, max_size: int=np.inf, processes: int=4, plotting: bool=False, verbose: bool=False, symmetric: bool=True, seed: int=0):
    if seed == -1:
        seed = random.randint(-10000000, 100000000)
    
    signature.columns = [0,1]
    if permutations < 1000 and not symmetric:
        if verbose:
            print('Low numer of permutations can lead to inaccurate p-value estimation. Symmetric Gamma distribution enabled to increase accuracy.')
        symmetric = True
    elif permutations < 500:
        if verbose:
            print('Low numer of permutations can lead to inaccurate p-value estimation. Consider increasing number of permutations.')
        symmetric = True

    signature = signature.sort_values(1, ascending=False).set_index(0)
    signature = signature[~signature.index.duplicated(keep='first')]
    abs_signature = signature
    abs_signature.iloc[:,0] = np.abs(abs_signature.iloc[:,0])

    signature_map = {}
    for i,h in enumerate(signature.index):
        signature_map[h] = i

    f_alpha_pos, f_beta_pos, f_pos_ratio, ks_pos, ks_neg = estimate_parameters(signature, abs_signature, signature_map, library, permutations=permutations, calibration_anchors=anchors, processes=processes, symmetric=symmetric, plotting=plotting, verbose=verbose, seed=seed)
    gsets = []
    
    params = []
    keys = list(library.keys())
    signature_genes = set(signature.index)
    for k in keys:
        stripped_set = strip_gene_set(signature, signature_genes, library[k])
        if len(stripped_set) >= min_size and len(stripped_set) <= max_size:
            gsets.append(k)
            params.append((signature, abs_signature, signature_map, stripped_set, f_alpha_pos, f_beta_pos, f_pos_ratio))
    
    with multiprocessing.Pool(processes) as pool:
        results = list(tqdm(pool.imap(probability_star, params), desc="Enrichment", total=len(params)))
    
    ess = []
    pvals = [1, 1]
    nes = []
    set_size = []
    legenes = []

    for res in results:
        gsize, es, ne, pval, legene = res
        nes.append(ne)
        ess.append(es)
        pvals.append(float(pval))
        set_size.append(gsize)
        legenes.append(legene)
    
    if not verbose:
        np.seterr(divide = 'ignore')
    
    fdr_values = multipletests(pvals, method="fdr_bh")[1]
    sidak_values = multipletests(pvals, method="sidak")[1]

    res =  pd.DataFrame([gsets, np.array(ess).astype("float"), np.array(nes).astype("float"), np.array(pvals).astype("float"), np.array(sidak_values).astype("float"), np.array(fdr_values).astype("float"), np.array(set_size).astype("int"), np.array(legenes)]).T
    res.columns = ["Term", "es", "nes", "pval", "sidak", "fdr","geneset_size", "leading_edge"]
    res = res.set_index("Term")

    if (ks_pos < 0.05 or ks_neg < 0.05) and verbose:
        print('Kolmogorov-Smirnov test failed. Gamma approximation deviates from permutation samples.\n'+"KS p-value (pos): "+str(ks_pos)+"\nKS p-value (neg): "+str(ks_neg))
    
    return res.sort_values("pval", key=abs, ascending=True)