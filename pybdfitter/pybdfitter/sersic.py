"""
Routines for dealing with Sersic profile
"""
from pylab import *
import numpy as np
from scipy.special import gamma
from .tools import nanomaggie2mag

__all__ = ["ip", "Sersic", "NSersic", "bn", "fn"]


# parameter index dictionary as used in bdfitter
ip = {
    'Ie': 0,
    'Re': 1,
    'n': 2,
    'q': 3,
    'c': 4,
    'x': 5,
    'y': 6,
    'pa': 7,
}

def bn(n):
    return n*exp(0.6950 - 0.1789/n)

def fn(n):
    return n * exp(bn(n)) * gamma(2.*n) / bn(n)**(2.*n)



class Sersic(object):
    """ 2D Sersic profile """
    param_names = ['Ie', 'Re', 'n', 'q', 'c', 'x', 'y', 'pa']

    def __init__(self, params):
        self.Ie = params[0]
        self.Re = params[1]
        self.n = params[2]
        self.q = params[3]
        self.c = params[4]
        self.x = params[5]
        self.y = params[6]
        self.pa = params[7]
        self.p = params

    # def __repr__(self):
    #     return '\n'.join([ "%2s = %g" % (n, self.__getattribute__(n)) \
    #             for n in Sersic.param_names ])

    @property
    def Ro(self):
        """ circularized radius """
        return self.Re*sqrt(self.q)

    @property
    def mu_e(self):
        """ surface brightness at r=Re (mag/arcsec^2) """
        return nanomaggie2mag(self.Ie/0.396**2)

    @property
    def bn(self):
        """ from Lima-Neto et al. 1999 """
        return self.n*exp(0.6950 - 0.1789/self.n)

    @property
    def mu_e_avg(self):
        """ mean surface brightness within Re mag/arcsec^2 """
        fn = self.n * exp(self.bn) * gamma(2.*self.n) / self.bn**(2.*self.n)
        return self.mu_e - 2.5*log10(fn)

    @property
    def total_flux_mag(self):
        """total flux in magnitude """
        return self.mu_e_avg - 2.5*log10(2*pi*(self.Re*0.396)**2*self.q)


class NSersic(object):
    """ N-Sersic profile"""
    def __init__(self, params):
        assert len(params) % 8 == 0, "params should be 8 X nprofiles"
        self.nprofiles = len(params) / 8
        self.profiles = []
        for i in range(self.nprofiles):
            self.profiles.append(Sersic(params[i*8:(i+1)*8]))

    def __getitem__(self, key):
        return self.profiles[key]

    def info(self):
        s = "NComponents = %4i" % (self.nprofiles)
        print s


def radius(param, x, y):
    iq = 1./param[0]
    c, x0, y0, phi = param[1:5]

    axtrans = np.abs((x-x0)*np.cos(phi) + (y-y0)*np.sin(phi))
    aytransq = np.abs(((y-y0) * np.cos(phi) - (x-x0) * np.sin(phi))*iq)
    if c > 1e-12:
        expo = c+2.0
        expo_inv = 1./expo
        r = (axtrans**expo + aytransq**expo)**expo_inv
    else:
        r = np.sqrt(axtrans**2 + aytransq**2)
    return r


def radius_fast(param, x, y):
    """ c = 0 """
    iq = 1./param[0]
    xx0 = x - param[2]
    yy0 = y - param[3]
    cp = np.cos(param[4])
    sp = np.sin(param[4])

    axtrans = xx0 * cp + yy0 * sp
    aytransq = (yy0*cp - xx0*sp)*iq
    r = np.sqrt(axtrans**2 + aytransq**2)
    return r


def sersic2d(x, y, p, cutoff=False):
    """
    get 1d or 2d Sersic model

    x, y : 1d or 2d array containing x, y positions
    p : Sersic profile parameter
    """

    Ie, reff, n, q, c, x0, y0, phi = p
    density = np.zeros_like(x)
    # non-sensical parameters
    if abs(Ie) < 1e-18 or reff < 1e-12:
        return density
    # get radius
    if abs(c) < 1e-12:
        r = radius_fast(p[3:], x, y)
    else:
        r = radius(p[3:], x, y)

    k = n*np.exp(0.6950 - 0.1789/n)
    ir0 = 1./reff

    # check if n is an easy number to work with
    if abs(n - 1.) < 1e-12:
        exponent = -k * (r*ir0 - 1.)
    elif abs(n - 4.) < 1e-12:
        exponent = -k*(np.sqrt(np.sqrt(r*ir0)) - 1.)
    elif abs(n - 2.) < 1e-12:
        exponent = -k*(np.sqrt(r*ir0) - 1.)
    else:
        exponent = -k * ((r*ir0)**(1./n) - 1.)
    density = Ie * np.exp(exponent)

    # cut the profile off
    # smoothly past a certain size
    # use the SDSS defined cutoffs
    if cutoff:
        if n < 1.05:
            fade = 3.*reff
            cut = 4.*reff
        else:
            fade = 7.*reff
            cut = 8.*reff
        denom = 1./(cut - fade)**2

        if len(x.shape) >= 1:
            truncate = np.where((r > fade) & (r < cut))
            cutoff = np.where(r > cut)
            if truncate:
                density[truncate] *= (1. - (r[truncate]-fade)**2*denom)**2
            if cutoff:
                density[cutoff] = 0.
        else:
            if r > cut:
                density = 0.
            if r > fade:
                density *= (1.-(r-fade)**2*denom)**2
    return density


def test_sersic2d():
    import pylab as pl
    x, y = np.meshgrid(np.arange(500), np.arange(500))
    param = [.5, 80., 4., 0.3, 0., 250., 250., 1.]
    s1 = sersic2d(x, y, param)
    param = [2., 20., 4., 1., 0., 250., 250., 1.]
    s2 = sersic2d(x, y, param)
    pl.imshow(np.log10(s1+s2))
    pl.colorbar()

    from pybdfitter.plot import get_profile
    pr = get_profile(s1, 1., 0, 250, 250, limit=200)
    pr2 = get_profile(s2, 1., 0, 250, 250, limit=200)
    pr_sum = get_profile(s1+s2, 1., 0, 250, 250, limit=200)

    pl.figure()
    pl.plot(pr.radius, pr.mnflux)
    pl.plot(pr2.radius, pr2.mnflux)
    pl.plot(pr_sum.radius, pr_sum.mnflux)
