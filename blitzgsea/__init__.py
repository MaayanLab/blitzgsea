import random
import numpy as np
import pandas as pd
from loess.loess_1d import loess_1d
from collections import Counter
from scipy import interpolate
from scipy.stats import norm
from scipy.stats import gamma
from scipy.stats import kstest
from matplotlib import pyplot as plt
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
import warnings

def strip_gene_set(signature, gene_set):
    signature_genes = set(signature.index)
    return [x for x in gene_set if x in signature_genes]

def enrichment_score(signature, signature_map, gene_set):
    gs = set(gene_set)
    hits = [signature_map[x] for x in gs if x in signature_map.keys()]
    hit_indicator = np.zeros(signature.shape[0])
    hit_indicator[hits] = 1
    no_hit_indicator = 1 - hit_indicator
    number_hits = len(hits)
    number_miss = signature.shape[0] - number_hits
    sum_hit_scores = np.sum(np.abs(signature.iloc[hits]))
    norm_hit =  float(1.0/sum_hit_scores)
    norm_no_hit = float(1.0/number_miss)
    running_sum = np.cumsum(hit_indicator * np.abs(signature.iloc[:,0]) * norm_hit - no_hit_indicator * norm_no_hit)
    nn = np.where(np.abs(running_sum)==np.max(np.abs(running_sum)))[0][0]
    es = running_sum[nn]
    return running_sum, es

def get_peak_size(signature, signature_map, size, permutations = 100):
    es = []
    for _ in range(permutations):
        rgenes = random.sample(list(signature.index), size)
        es.append(enrichment_score(signature, signature_map, rgenes)[1])
    return es

def loess_interpolation(x, y):
    yl = np.array(y)
    xout, yout, wout = loess_1d(x, yl)
    return interpolate.interp1d(xout, yout)

def estimate_parameters(signature, signature_map, library, permutations: int=1000, symmetric: bool=False, plotting: bool=False):
    ll = []
    for key in library.keys():
        ll.append(len(library[key]))
    cc = Counter(ll)
    set_sizes = pd.DataFrame(list(cc.items()),columns = ['set_size','count']).sort_values("set_size")
    set_sizes["cumsum"] = np.cumsum(set_sizes.iloc[:,1])
    
    x = [1,2,3,4,5,8,10,12,14,18,22,26,30,35,60,100,200,500, set_sizes.iloc[:,0].max()]
         
    alpha_pos = []
    beta_pos = []
    loc_pos = []

    alpha_neg = []
    beta_neg = []
    loc_neg = []

    pos_ratio = []

    ks_pos = []
    ks_neg = []
    
    pbar = tqdm(x)
    for i in pbar:
        pbar.set_description("Parameter Calibration %s" % i)
        es = np.array(get_peak_size(signature, signature_map, i, permutations=permutations))
        pos = [x for x in es if x > 0]
        neg = [x for x in es if x < 0]
        if (len(neg) < 250 or len(pos) < 250) and not symmetric:
            print('Low numer of permutations can lead to inaccurate p-value estimation. Symmetric Gamma distribution enabled to increase accuracy.')
            symmetric = True
        if symmetric:
            aes = np.abs(es)
            fit_alpha, fit_loc, fit_beta = gamma.fit(aes)
            ks_pos.append(kstest(aes, 'gamma', args=(fit_alpha, fit_loc, fit_beta))[1])
            ks_neg.append(kstest(aes, 'gamma', args=(fit_alpha, fit_loc, fit_beta))[1])
            alpha_pos.append(fit_alpha)
            beta_pos.append(fit_beta)
            loc_pos.append(fit_loc)
            alpha_neg.append(fit_alpha)
            beta_neg.append(fit_beta)
            loc_neg.append(fit_loc)
        else:
            fit_alpha, fit_loc, fit_beta = gamma.fit(pos)
            ks_pos.append(kstest(pos, 'gamma', args=(fit_alpha, fit_loc, fit_beta))[1])
            alpha_pos.append(fit_alpha)
            beta_pos.append(fit_beta)
            loc_pos.append(fit_loc)
            
            fit_alpha, fit_loc, fit_beta = gamma.fit(-np.array(neg))
            ks_neg.append(kstest(-np.array(neg), 'gamma', args=(fit_alpha, fit_loc, fit_beta))[1])
            alpha_neg.append(fit_alpha)
            beta_neg.append(fit_beta)
            loc_neg.append(fit_loc)
        pos_ratio.append(len(pos)/(len(pos)+len(neg)))
    pbar.close()

    if np.max(pos_ratio) > 1.5:
        print('Significant unbalance between positive and negative enrichment scores detected. Signature values are not centered close to 0.')

    x = np.array(x, dtype=float)
    
    f_alpha_pos = loess_interpolation(x, alpha_pos)
    f_beta_pos = loess_interpolation(x, beta_pos)
    f_loc_pos = loess_interpolation(x, loc_pos)
    
    f_alpha_neg = loess_interpolation(x, alpha_neg)
    f_beta_neg = loess_interpolation(x, beta_neg)
    f_loc_neg = loess_interpolation(x, loc_neg)

    f_pos_ratio = loess_interpolation(x, pos_ratio)
    
    if plotting:
        xx = np.linspace(min(x), max(x), 1000)
        
        plt.figure(0)
        yy = f_loc_pos(xx)
        plt.plot(xx, yy, '--', lw=3)
        plt.plot(x, loc_pos, 'ko')
        
        yy = f_loc_neg(xx)
        plt.plot(xx, yy, '--', lw=3, c="orange")
        plt.plot(x, loc_neg, 'o', c="coral")

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
        
    return f_alpha_pos, f_beta_pos, f_loc_pos, f_alpha_neg, f_beta_neg, f_loc_neg, f_pos_ratio, np.mean(ks_pos), np.mean(ks_neg)

def gsea(signature, library, permutations: int=100, plotting: bool=False, verbose: bool=False, symmetric: bool=False):
    if permutations < 1000 and not symmetric:
        print('Low numer of permutations can lead to inaccurate p-value estimation. Symmetric Gamma distribution enabled to increase accuracy.')
        symmetric = True
    elif permutations < 500:
        print('Low numer of permutations can lead to inaccurate p-value estimation.')
        symmetric = True

    signature = signature.sort_values(1, ascending=False).set_index(0)
    signature_map = {}
    for i,h in enumerate(signature.index):
        signature_map[h] = i

    f_alpha_pos, f_beta_pos, f_loc_pos, f_alpha_neg, f_beta_neg, f_loc_neg, f_pos_ratio, ks_pos, ks_neg = estimate_parameters(signature, signature_map, library, permutations=permutations, symmetric=symmetric, plotting=plotting)
    gsets = []
    ess = []
    pvals = []
    nes = []
    set_size = []

    lib_keys = list(library.keys())
    pbar = tqdm(range(len(lib_keys)))
    
    for gene_set_key_i in pbar:
        gene_set_key = lib_keys[gene_set_key_i]
        pbar.set_description("GSEA %s" % gene_set_key_i)
        gene_set = strip_gene_set(signature, library[gene_set_key])
        gsize = len(gene_set)
        if gsize > 0:
            set_size.append(gsize)
            rs, es = enrichment_score(signature, signature_map, gene_set)

            pos_alpha = f_alpha_pos(gsize)
            pos_beta = f_beta_pos(gsize)
            pos_loc = f_loc_pos(gsize)
            neg_alpha = f_alpha_neg(gsize)
            neg_beta = f_beta_neg(gsize)
            neg_loc = f_loc_neg(gsize)

            pos_ratio = f_pos_ratio(gsize)

            gsets.append(gene_set_key)
            ess.append(es)

            if es > 0:
                rv = gamma(pos_alpha, scale=pos_beta, loc=pos_loc)
                prob = rv.cdf(es)
                prob_two_tailed = np.min([0.5,(1-np.min([(1-pos_ratio)+prob*pos_ratio,1]))])
                if prob_two_tailed == 1:
                    nes.append(0)
                else:
                    nes.append(norm.ppf(1-np.min([1,prob_two_tailed])))
                pvals.append(2*prob_two_tailed)
            else:
                rv = gamma(neg_alpha, scale=neg_beta, loc=neg_loc)
                prob = rv.cdf(-es)
                prob_two_tailed = np.min([0.5,(1-np.min([prob*(1-pos_ratio)+pos_ratio,1]))])
                nes.append(norm.ppf(np.min([1,prob_two_tailed])))
                pvals.append(2*prob_two_tailed)
    pbar.close()

    fdr_values = multipletests(pvals, method="fdr_bh")[1]
    sidak_values = multipletests(pvals, method="sidak")[1]

    res =  pd.DataFrame([gsets, np.array(ess).astype("float"), np.array(nes).astype("float"), np.array(pvals).astype("float"), np.array(sidak_values).astype("float"), np.array(fdr_values).astype("float"), np.array(set_size).astype("int")]).T
    res.columns = ["Term", "es", "zscore", "pval", "sidak", "fdr","geneset_size"]
    res = res.set_index("Term")

    if ks_pos < 0.05 or ks_neg < 0.05:
        print('Kolmogorov-Smirnov test failed. Gamma approximation deviates from permutation samples.\n'+"KS p-value (pos): "+str(ks_pos)+"\nKS p-value (neg): "+str(ks_neg))
    
    return res.sort_values("pval", key=abs, ascending=True)