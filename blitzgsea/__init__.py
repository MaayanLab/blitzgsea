import random
import numpy as np
import pandas as pd
from loess.loess_1d import loess_1d
from statsmodels.nonparametric.smoothers_lowess import lowess
from collections import Counter
from scipy import interpolate

from matplotlib import pyplot as plt
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
import multiprocessing

from mpmath import mp
from mpmath import mpf
from mpsci.distributions.normal import invcdf
from mpsci.distributions.gamma import cdf as gammacdf

from scipy.stats import gamma
from scipy.stats import kstest

import blitzgsea.signature_similarity
from blitzgsea.signature_similarity import create_pdf, map_density_range, kl_divergence, best_kl_fit
import blitzgsea.enrichr
import blitzgsea.plot
import blitzgsea.shuffle

from importlib import reload
reload(blitzgsea.signature_similarity)
reload(blitzgsea.enrichr)
reload(blitzgsea.plot)
reload(blitzgsea.shuffle)

mp.dps = 1000
mp.prec = 1000
global pdf_cache
pdf_cache = {}

def strip_gene_set(signature, signature_genes, gene_set):
    return [x for x in gene_set if x in signature_genes]

def enrichment_score(abs_signature, signature_map, gene_set):
    hits = [signature_map[x] for x in gene_set if x in signature_map.keys()]
    hit_indicator = np.zeros(len(abs_signature))
    hit_indicator[hits] = 1
    no_hit_indicator = 1 - hit_indicator
    number_hits = len(hits)
    number_miss = len(abs_signature) - number_hits
    sum_hit_scores = np.sum(abs_signature[hits])
    if sum_hit_scores == 0:
        return 0, 0
    norm_hit = float(1.0/sum_hit_scores)
    norm_no_hit = float(1.0/number_miss)
    running_sum = np.cumsum(hit_indicator * abs_signature * norm_hit - no_hit_indicator * norm_no_hit)
    nn = np.argmax(np.abs(running_sum))
    es = running_sum[nn]
    return running_sum, es

def enrichment_score_null(abs_signature, hit_indicator, number_hits):
    np.random.shuffle(hit_indicator)
    hits = np.where(hit_indicator == 1)[0]
    no_hit_indicator = 1 - hit_indicator
    number_miss = len(abs_signature) - number_hits
    sum_hit_scores = np.sum(abs_signature[hits])
    if sum_hit_scores == 0:
        return 0
    norm_hit = float(1.0/sum_hit_scores)
    norm_no_hit = float(1.0/number_miss)
    running_sum = np.cumsum(hit_indicator * abs_signature * norm_hit - no_hit_indicator * norm_no_hit)
    peak = np.abs(running_sum).argmax()
    es = running_sum[peak]
    return es

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
        es.append(enrichment_score(abs_signature, signature_map, rgenes)[1])
    return es

def get_peak_size_adv(abs_signature, size, permutations, seed):
    random.seed(seed)
    es = []
    number_hits = size
    hit_indicator = np.zeros(len(abs_signature))
    hit_indicator[0:number_hits] = 1
    for i in range(permutations):
        es.append(enrichment_score_null(abs_signature, hit_indicator, number_hits))
    return es

def loess_interpolation(x, y, frac=0.5):
    yl = np.array(y)
    #xout, yout, wout = loess_1d(x, yl, frac=frac)
    #return interpolate.interp1d(xout, yout)
    yout = lowess(y, x, frac=frac)[:, 1]
    return interpolate.interp1d(x, yout)

def estimate_parameters(signature, abs_signature, signature_map, library, permutations: int=2000, max_size=4000, symmetric: bool=False, calibration_anchors: int=20, plotting: bool=False, processes=4, verbose=False, progress=False, seed: int=0):
    ll = []
    for key in library.keys():
        ll.append(len(library[key]))
    cc = Counter(ll)
    set_sizes = pd.DataFrame(list(cc.items()),columns = ['set_size','count']).sort_values("set_size")
    set_sizes["cumsum"] = np.cumsum(set_sizes.iloc[:,1])
    
    ll = [len(library[l]) for l in library]
    nn = np.percentile(ll, q=np.linspace(2, 100, calibration_anchors))
    anchor_set_sizes = sorted(list(set(np.append([1,4,6, np.min([max_size, np.max(ll)]), np.min([max_size, int(signature.shape[0]/2)]), np.min([max_size, signature.shape[0]-1])], nn).astype("int"))))
    anchor_set_sizes = [x for x in anchor_set_sizes if x < signature.shape[0]]

    with multiprocessing.Pool(processes) as pool:
        args = [(signature, abs_signature, signature_map, xx, permutations, symmetric, seed+xx) for xx in anchor_set_sizes]
        results = list(tqdm(pool.imap(estimate_anchor_star, args), desc="Calibration", total=len(args), disable=not progress))

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

    anchor_set_sizes = np.array(anchor_set_sizes, dtype=float)
    
    f_alpha_pos = loess_interpolation(anchor_set_sizes, alpha_pos)
    f_beta_pos = loess_interpolation(anchor_set_sizes, beta_pos, frac=0.2)
    
    f_alpha_neg = loess_interpolation(anchor_set_sizes, alpha_neg)
    f_beta_neg = loess_interpolation(anchor_set_sizes, beta_neg, frac=0.2)

    # fix issue with numeric instability
    pos_ratio = pos_ratio - np.abs(0.0001*np.random.randn(len(pos_ratio)))
    f_pos_ratio = loess_interpolation(anchor_set_sizes, pos_ratio)
    
    if plotting:
        xx = np.linspace(min(anchor_set_sizes), max(anchor_set_sizes), 1000)
        
        plt.figure(0)
        
        yy = f_alpha_pos(xx)
        plt.figure(1)
        plt.plot(xx, yy, '--', lw=3)
        plt.plot(anchor_set_sizes, alpha_pos, 'ko')
        
        yy = f_alpha_neg(xx)
        plt.figure(1)
        plt.plot(xx, yy, '--', lw=3, c="orange")
        plt.plot(anchor_set_sizes, alpha_neg, 'o', c="coral")

        yy = f_beta_pos(xx)
        plt.figure(1)
        plt.plot(xx, yy, '--', lw=3)
        plt.plot(anchor_set_sizes, beta_pos, 'ko')
        
        yy = f_beta_neg(xx)
        plt.figure(1)
        plt.plot(xx, yy, '--', lw=3, c="orange")
        plt.plot(anchor_set_sizes, beta_neg, 'o', c="coral")
        
        yy = f_pos_ratio(xx)
        plt.figure(2)
        plt.plot(xx, yy, lw=3)
        plt.plot(anchor_set_sizes, pos_ratio, 'o', c="black")
        
    return f_alpha_pos, f_beta_pos, f_pos_ratio, np.mean(ks_pos), np.mean(ks_neg)

def estimate_anchor_star(args):
    return estimate_anchor(*args)

def estimate_anchor(signature, abs_signature, signature_map, set_size, permutations, symmetric, seed):
    es = np.array(get_peak_size_adv(abs_signature, set_size, permutations, seed))
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

def gsea(signature, library, permutations: int=2000, anchors: int=20, min_size: int=5, max_size: int=4000, processes: int=4, plotting: bool=False, verbose: bool=False, progress: bool=False, symmetric: bool=True, signature_cache: bool=True, kl_threshold: float=0.3, kl_bins: int=200, shared_null: bool=False, seed: int=0, add_noise: bool=False):
    """
    Perform Gene Set Enrichment Analysis (GSEA) on the given signature and library.

    Parameters:
    signature (array-like): The gene expression signature to analyze.
    library (array-like): The gene set library to use for enrichment analysis.
    permutations (int, optional): Number of randomized permutations to estimate ES distributions. Default is 2000.
    anchors (int, optional): Number of gene set size distributions calculated. Remaining are interpolated. Default is 20.
    min_size (int, optional): Minimum number of genes in geneset. Default is 5.
    max_size (int, optional): Maximal number of genes in gene set. Default is 4000.
    processes (int, optional): Number of parallel threads. Not much gain after 4 threads. Default is 4.
    symmetric (bool, optional): Use same distribution parameters for negative and positive ES. If False estimate them separately. Default is True.
    signature_cache (bool, optional): Use precomputed anchor parameters. Default is True.
    shared_null (bool, optional): Use same null for signatures if a compatible model already exists (uses KL-divergence test). Default is False.
    kl_threshold (float, optional): Controls how similar signature value distributions have to be for reuse. Default is 0.3.
    kl_bins (int, optional): Number of bins in PDF representation of distributions for KL test. Default is 200.
    plotting (bool, optional): Plot estimated anchor parameters. Default is False.
    verbose (bool, optional): Toggle additional output. Default is False.
    progress (bool, optional): Toggle progress bar. Default is False.
    seed (int, optional): Random seed. Same seed will result in identical results. If seed equals -1, generate a random seed. Default is 0.
    add_noise (bool, optional): Add small random noise to signature. The noise is a fraction of the expression values. Default is False.

    Returns:
    array-like: Enrichment scores for each gene set in the library.
    """
    if seed == -1:
        seed = random.randint(-10000000, 100000000)

    signature.columns = ["i","v"]
    if permutations < 1000 and not symmetric:
        if verbose:
            print('Low numer of permutations can lead to inaccurate p-value estimation. Symmetric Gamma distribution enabled to increase accuracy.')
        symmetric = True
    elif permutations < 500:
        if verbose:
            print('Low numer of permutations can lead to inaccurate p-value estimation. Consider increasing number of permutations.')
        symmetric = True

    random.seed(seed)
    sig_hash = hash(signature.to_string())

    # optionally noise can be added as a fraction of the expression values
    if add_noise:
        signature.iloc[:,1] = signature.iloc[:,1] + np.random.normal(signature.shape[0])/(np.mean(np.abs(signature.iloc[:,1]))*100000)
    signature = signature.sort_values("v", ascending=False).set_index("i")
    signature = signature[~signature.index.duplicated(keep='first')]
    
    abs_signature = np.array(np.abs(signature.iloc[:,0]))

    signature_map = {}
    for i,h in enumerate(signature.index):
        signature_map[h] = i

    if shared_null and len(pdf_cache) > 0:
        kld, sig_hash_temp = best_kl_fit(np.array(signature["v"]), pdf_cache, bins=kl_bins)
        if kld < kl_threshold:
            sig_hash = sig_hash_temp
            if verbose:
                print(f"Found compatible null model. Best KL-divergence: {kld}")
        elif verbose:
            print(f"No compatible null model. Best KL-divergence: {kld} > kl_threshold: {kl_threshold}")

    if sig_hash in pdf_cache.keys() and signature_cache:
        if verbose:
            print("Use cached anchor parameters")
        f_alpha_pos, f_beta_pos, f_pos_ratio, ks_pos, ks_neg = pdf_cache[sig_hash]["model"]
    else:
        f_alpha_pos, f_beta_pos, f_pos_ratio, ks_pos, ks_neg = estimate_parameters(signature, abs_signature, signature_map, library, permutations=permutations, calibration_anchors=anchors, processes=processes, symmetric=symmetric, plotting=plotting, verbose=verbose, seed=seed, progress=progress, max_size=max_size)
        pdf_cache[sig_hash] = {}
        xv, pdf = create_pdf(np.array(signature["v"]), kl_bins)
        pdf_cache[sig_hash]["xvalues"] = xv
        pdf_cache[sig_hash]["pdf"] = pdf
        pdf_cache[sig_hash]["model"] = (f_alpha_pos, f_beta_pos, f_pos_ratio, ks_pos, ks_neg)
    
    gsets = []

    keys = list(library.keys())
    signature_genes = set(signature.index)

    ess = []
    pvals = []
    ness = []
    set_size = []
    legeness = []
    
    for k in tqdm(keys, desc="Enrichment ", disable=not verbose):
        stripped_set = strip_gene_set(signature, signature_genes, library[k])
        if len(stripped_set) >= min_size and len(stripped_set) <= max_size:
            gsets.append(k)
            gsize = len(stripped_set)
            rs, es = enrichment_score(abs_signature, signature_map, stripped_set)
            legenes = get_leading_edge(rs, signature, stripped_set, signature_map)

            pos_alpha = f_alpha_pos(gsize)
            pos_beta = f_beta_pos(gsize)
            pos_ratio = f_pos_ratio(gsize)

            mp.dps = 50
            mp.prec = 50

            if es > 0:
                prob = gamma.cdf(es, float(pos_alpha), scale=float(pos_beta))
                if prob > 0.99999999:
                    mp.dps = 500
                    mp.prec = 500
                    prob = gammacdf(es, float(pos_alpha), float(pos_beta))
                #prob_two_tailed = np.min([0.5,(1-np.min([(1-pos_ratio)+prob*pos_ratio,1]))])
                prob_two_tailed = np.min([0.5,(1-np.min([prob*pos_ratio+1-pos_ratio,1]))])
                if prob_two_tailed == 1:
                    nes = 0
                else:
                    nes = invcdf(mpf(1)-mpf(np.min([1,prob_two_tailed])))
                pval = 2*prob_two_tailed
            else:
                prob = gamma.cdf(-es, float(pos_alpha), scale=float(pos_beta))
                if prob > 0.99999999:
                    mp.dps = 500
                    mp.prec = 500
                    prob = gammacdf(-es, float(pos_alpha), float(pos_beta))
                # prob_two_tailed = np.min([0.5,(1-np.min([prob*(1-pos_ratio)+pos_ratio,1]))])
                prob_two_tailed = np.min([0.5,(1-np.min([(((prob)-(prob*pos_ratio))+pos_ratio),1]))])
                nes = invcdf(mpf(np.min([1,prob_two_tailed])))
                pval = 2*prob_two_tailed
            
            mp.dps = 50
            mp.prec = 50

            ness.append(float(nes))
            ess.append(float(es))
            pvals.append(float(pval))
            set_size.append(gsize)
            legeness.append(legenes)

    if not verbose:
        np.seterr(divide = 'ignore')
    
    fdr_values = multipletests(pvals, method="fdr_bh")[1]
    sidak_values = multipletests(pvals, method="sidak")[1]

    res =  pd.DataFrame([gsets, np.array(ess), np.array(ness), np.array(pvals), np.array(sidak_values), np.array(fdr_values), np.array(set_size), np.array(legeness)]).T
    res.columns = ["Term", "es", "nes", "pval", "sidak", "fdr","geneset_size", "leading_edge"]
    res["Term"] = res['Term'].astype("str")
    res["es"] = res['es'].astype("float")
    res["nes"] = res['nes'].astype("float")
    res["pval"] = res['pval'].astype("float")
    res["sidak"] = res['sidak'].astype("float")
    res["fdr"] = res['fdr'].astype("float")
    res["geneset_size"] = res['geneset_size'].astype("int")
    res["leading_edge"] = res['leading_edge'].astype("str")
    res = res.set_index("Term")

    if (ks_pos < 0.05 or ks_neg < 0.05) and verbose:
        print('Kolmogorov-Smirnov test failed. Gamma approximation deviates from permutation samples.\n'+"KS p-value (pos): "+str(ks_pos)+"\nKS p-value (neg): "+str(ks_neg))
    
    return res.sort_values("pval", key=abs, ascending=True)
