import random
import numpy as np
import pandas as pd
import loess
from collections import Counter
from scipy import interpolate
from scipy.stats import norm
from matplotlib import pyplot as plt


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
    
    x = [4,5,8,10,12,14,18,22,26,30,35]
         
    means_pos = []
    sd_pos = []
    
    means_neg = []
    sd_neg = []
    
    pos_ratio = []
    
    for i in x:
        print(i)
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

