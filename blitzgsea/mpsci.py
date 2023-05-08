# this is directly copied from the fantastic project https://github.com/WarrenWeckesser/mpsci

from mpmath import mp

def gammacdf(x, k, theta):
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

def invcdf(p, mu=0, sigma=1):
    """
    Normal distribution inverse CDF.
    This function is also known as the quantile function or the percent
    point function.
    """
    p = min(max(p, 0), 1)
    with mp.extradps(mp.dps):
        mu = mp.mpf(mu)
        sigma = mp.mpf(sigma)

        a = mp.erfinv(2*p - 1)
        x = mp.sqrt(2)*sigma*a + mu
        return x
