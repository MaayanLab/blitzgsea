import numpy as np
import pandas as pd

from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d

def create_pdf(data, bins=200):
    x_vals = np.linspace(data.min(), data.max(), num=bins)
    kde = gaussian_kde(data)
    return x_vals, kde(x_vals)

def map_density_range(new_range, old_range, old_pdf):
    pdf_interp = interp1d(old_range, old_pdf, kind='linear', fill_value=0.0, bounds_error=False)
    return pdf_interp(new_range)

def kl_divergence(p, q):
    is_finite = (p != 0) & (q != 0)
    p, q = p[is_finite], q[is_finite]
    return np.sum(p * (np.log2(p) - np.log2(q)))

def best_kl_fit(signature, pdfs, bins=200):
    xv, pdf = create_pdf(signature, bins)
    res = []
    pdf_keys = list(pdfs.keys())
    for k in pdf_keys:
        ipdf = map_density_range(xv, pdfs[k]["xvalues"], pdfs[k]["pdf"])
        k = (kl_divergence(pdf, ipdf)+kl_divergence(ipdf,pdf))/2
        res.append(k)
    min_pos = np.argmin(res)
    return res[min_pos], pdf_keys[min_pos]
