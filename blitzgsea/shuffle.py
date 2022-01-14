import numpy as np
from statsmodels.stats.multitest import multipletests
import pandas as pd

def gsea(exprs, library, groups, permutations=1000, seed=1):

    pos = 0
    neg = 1

    rs = np.random.RandomState(seed)
    expr_mat = exprs.T
    perm_cor_tensor = np.tile(expr_mat, (permutations,1,1))

    perm_cor_tensor.shape
    for arr in perm_cor_tensor[:-1]: rs.shuffle(arr)

    groups = np.array(groups)
    pos = groups == pos
    neg = groups == neg
    n_pos = np.sum(pos)
    n_neg = np.sum(neg)
    pos_cor_mean = perm_cor_tensor[:,pos,:].mean(axis=1)
    neg_cor_mean = perm_cor_tensor[:,neg,:].mean(axis=1)
    pos_cor_std = perm_cor_tensor[:,pos,:].std(axis=1, ddof=1)
    neg_cor_std = perm_cor_tensor[:,neg,:].std(axis=1, ddof=1)

    denom = np.sqrt((pos_cor_std**2)/n_pos  + (neg_cor_std**2)/n_neg)
    cor_mat = (pos_cor_mean - neg_cor_mean)/ denom

    cor_mat_ind = cor_mat.argsort()
    gene_mat = cor_mat_ind[:, ::-1]
    cor_mat = cor_mat[:, ::-1]

    cor_mat = cor_mat.T
    keys = np.array(list(library.keys()))

    genes = np.array(exprs.index)
    genes_ind = gene_mat

    tag_indicator = np.vstack([np.in1d(genes, library[key], assume_unique=True) for key in keys])
    tag_indicator = tag_indicator.astype(int)
    perm_tag_tensor = np.stack([tag.take(genes_ind).T for tag in tag_indicator], axis=0)
    
    no_tag_tensor = 1 - perm_tag_tensor
    rank_alpha = np.abs(perm_tag_tensor*cor_mat[np.newaxis,:,:])

    axis=1
    P_GW_denominator = np.sum(rank_alpha, axis=axis, keepdims=True)
    P_NG_denominator = np.sum(no_tag_tensor, axis=axis, keepdims=True)
    REStensor = np.cumsum(rank_alpha / P_GW_denominator - no_tag_tensor / P_NG_denominator, axis=axis)

    esmax, esmin = REStensor.max(axis=axis), REStensor.min(axis=axis)
    esmatrix = np.where(np.abs(esmax) > np.abs(esmin), esmax, esmin)

    es, esnull, RES = esmatrix[:,-1], esmatrix[:,:-1], REStensor[:,:,-1]

    condlist = [ es < 0, es >=0]
    choicelist = [(esnull < es.reshape(len(es),1)).sum(axis=1)/ (esnull < 0).sum(axis=1),
                (esnull >= es.reshape(len(es),1)).sum(axis=1)/ (esnull >= 0).sum(axis=1)]
    
    pvals = np.select(condlist, choicelist)
    fdr_values = multipletests(pvals, method="fdr_bh")[1]
    sidak_values = multipletests(pvals, method="sidak")[1]

    res =  pd.DataFrame([keys, np.array(es).astype("float"), np.array(pvals).astype("float"), np.array(fdr_values).astype("float"), np.array(sidak_values).astype("float")]).T
    res.columns = ["Term", "es", "pval", "fdr", "sidak"]
    res = res.set_index("Term")

    return res.sort_values("pval")
