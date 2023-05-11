# this is directly copied from the fantastic project https://github.com/WarrenWeckesser/mpsci

from mpmath import mp
import math
from scipy.stats import norm, gamma


def gammacdf_old(x, k, theta):
    """
    Gamma distribution cumulative distribution function.
    k is the shape parameter
    theta is the scale parameter (reciprocal of the rate parameter)
    Unlike scipy, a location parameter is not included.
    """
    with mp.extradps(mp.dps):
        x = mp.mpf(x)
        if x < 0:
            return mp.zero
        return mp.gammainc(k, 0, x/theta, regularized=True)

def gammacdf(x, k, theta):
    log_sf_value = gamma.logsf(x, k, scale=theta)
    sf_value = mp.exp(log_sf_value)
    cdf_value = 1 - sf_value
    return cdf_value

def invcdf_old(p, mu=0, sigma=1):
    """
    Normal distribution inverse CDF.
    This function is also known as the quantile function or the percent
    point function.
    """
    if math.isnan(p):
        p = 1
    p = min(max(p, 0), 1)
    with mp.extradps(mp.dps):
        mu = mp.mpf(mu)
        sigma = mp.mpf(sigma)
        try:
            a = mp.erfinv(2*p - 1)
            x = mp.sqrt(2)*sigma*a + mu
        except Exception:
            print("The problem value is: ", p)
            quit()
        return x

def invcdf(p, mu=0, sigma=1):
    p = float(p)
    if math.isnan(p):
        p = 1
    p = min(max(p, 0), 1)
    n = norm.isf(p)
    return n