import random
import numpy as np
import pandas as pd
from loess.loess_1d import loess_1d
from collections import Counter
from scipy import interpolate
from scipy.stats import norm
from matplotlib import pyplot as plt
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests


def strip_gene_set(signature, gene_set):
    signature_genes = set(signature.iloc[:,0])
    return [x for x in gene_set if x in signature_genes]
    
def enrichment_score(signature, gene_set):
    rank_vector = signature.sort_values(1, ascending=False).set_index(0)
    gs = set(gene_set)
    hits = [i for i,x in enumerate(rank_vector.index) if x in gs]
    hit_indicator = np.zeros(rank_vector.shape[0])
    hit_indicator[hits] = 1
    no_hit_indicator = 1 - hit_indicator
    number_hits = len(hits)
    number_miss = rank_vector.shape[0] - number_hits
    sum_hit_scores = np.sum(np.abs(rank_vector.iloc[hits]))
    norm_hit =  1.0/sum_hit_scores
    norm_no_hit = 1.0/number_miss
    running_sum = np.cumsum(hit_indicator * np.abs(rank_vector[1]) * float(norm_hit) - no_hit_indicator * float(norm_no_hit))
    nn = np.where(np.abs(running_sum)==np.max(np.abs(running_sum)))[0][0]
    es = running_sum[nn]
    return running_sum, es

def get_peak_size(signature, size, permutations = 100):
    es = []
    for _ in range(permutations):
        rgenes = random.sample(list(signature.iloc[:,0]), size)
        es.append(enrichment_score(signature, rgenes)[1])
    return es

def loess_interpolation(x, y):
    yl = np.array(y)
    xout, yout, wout = loess_1d(x, yl)
    return interpolate.interp1d(xout, yout)

def estimate_parameters(signature, library, permutations: int=1000, plotting: bool=False):
    ll = []
    for key in library.keys():
        ll.append(len(library[key]))
    cc = Counter(ll)
    set_sizes = pd.DataFrame(list(cc.items()),columns = ['set_size','count']).sort_values("set_size")
    set_sizes["cumsum"] = np.cumsum(set_sizes.iloc[:,1])
    
    x = [1,2,3,4,5,8,10,12,14,18,22,26,30,35,60,100,200,500, set_sizes.iloc[:,0].max()]
         
    means_pos = []
    sd_pos = []
    means_neg = []
    sd_neg = []
    pos_ratio = []
    
    pbar = tqdm(x)
    for i in pbar:
        pbar.set_description("Parameter Calibration %s" % i)
        es = np.array(get_peak_size(signature, i, permutations=permutations))
        pos = [x for x in es if x > 0]
        param = norm.fit(pos)
        means_pos.append(param[0])
        sd_pos.append(param[1])
        neg = [x for x in es if x < 0]
        param = norm.fit(neg)
        means_neg.append(param[0])
        sd_neg.append(param[1])
        pos_ratio.append(len(pos)/(len(pos)+len(neg)))
    pbar.close()

    x = np.array(x, dtype=float)
    
    f_mean_pos = loess_interpolation(x, means_pos)
    f_sd_pos = loess_interpolation(x, sd_pos)
    
    f_mean_neg = loess_interpolation(x, means_neg)
    f_sd_neg = loess_interpolation(x, sd_neg)
    
    f_pos_ratio = loess_interpolation(x, pos_ratio)
    
    if plotting:
        xx = np.linspace(min(x), max(x), 1000)
        
        plt.figure(0)
        yy = f_mean_pos(xx)
        plt.plot(xx, yy, '--', lw=3)
        plt.plot(x, means_pos, 'ko')
        
        yy = f_mean_neg(xx)
        plt.plot(xx, yy, '--', lw=3, c="orange")
        plt.plot(x, means_neg, 'o', c="coral")

        yy = f_sd_pos(xx)
        plt.figure(1)
        plt.plot(xx, yy, '--', lw=3)
        plt.plot(x, sd_pos, 'ko')
        
        yy = f_sd_neg(xx)
        plt.figure(1)
        plt.plot(xx, yy, '--', lw=3, c="orange")
        plt.plot(x, sd_neg, 'o', c="coral")
        
        yy = f_pos_ratio(xx)
        plt.figure(2)
        plt.plot(xx, yy, lw=3)
        plt.plot(x, pos_ratio, 'o', c="black")
        
    return f_mean_pos, f_sd_pos, f_mean_neg, f_sd_neg, f_pos_ratio

def gsea(signature, library, permutations: int=100, plotting: bool=False, verbose: bool=False):
    f_mean_pos, f_sd_pos, f_mean_neg, f_sd_neg, f_pos_ratio = estimate_parameters(signature, library, permutations=permutations, plotting=plotting)
    gsets = []
    ess = []
    pvals = []
    nes = []

    lib_keys = list(library.keys())
    pbar = tqdm(range(len(lib_keys)))
    
    for gene_set_key_i in pbar:
        gene_set_key = lib_keys[gene_set_key_i]
        pbar.set_description("GSEA %s" % gene_set_key_i)
        gene_set = strip_gene_set(signature, library[gene_set_key])
        gsize = len(gene_set)
        if gsize > 0:
            rs, es = enrichment_score(signature, gene_set)

            pos_mean = f_mean_pos(gsize)
            pos_sd = f_sd_pos(gsize)
            neg_mean = f_mean_neg(gsize)
            neg_sd = f_sd_neg(gsize)

            pos_ratio = f_pos_ratio(gsize)

            gsets.append(gene_set_key)
            ess.append(es)

            if es > 0:
                rv = norm(loc=pos_mean, scale=pos_sd)
                prob = rv.cdf(es)
                prob_two_tailed = np.min([2*(1-((1-pos_ratio)+prob*pos_ratio)), 1])
                if prob_two_tailed == 1:
                    nes.append(0)
                else:
                    nes.append(norm.ppf(1-prob_two_tailed))
                pvals.append(prob_two_tailed)
            else:
                rv = norm(loc=neg_mean, scale=neg_sd)
                prob = 1-rv.cdf(es)
                prob_two_tailed = np.min([2*(1-(prob*(1-pos_ratio)+pos_ratio)),1])
                nes.append(norm.ppf(prob_two_tailed))
                pvals.append(prob_two_tailed)
    pbar.close()

    fdr_values = multipletests(pvals, method="fdr_bh")[1]
    sidak_values = multipletests(pvals, method="sidak")[1]

    res =  pd.DataFrame([gsets, ess, nes, pvals, sidak_values, fdr_values]).T
    res.columns = ["gene_set", "es", "zscore", "pval", "sidak", "fdr"]
    
    return res.sort_values("pval", key=abs, ascending=True)