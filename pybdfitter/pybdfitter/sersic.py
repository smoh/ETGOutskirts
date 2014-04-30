"""
Routines for dealing with Sersic profile
"""
from pylab import *
from scipy.special import gamma

__all__ = ["nanomaggie2mag", "ip", "Sersic"]

def nanomaggie2mag(nnmg):
    """ nnmg : flux in nanomaggie """
    return 22.5 - 2.5*log10(nnmg)

def mag2nanomaggie(mag):
    return 10**((22.5 - mag)/2.5)

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

    def __repr__(self):
        return '\n'.join([ "%2s = %g" % (n, self.__getattribute__(n)) \
                for n in Sersic.param_names ])
    
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
