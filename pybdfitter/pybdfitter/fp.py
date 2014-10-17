import numpy as np
from astropy.table import Table
from astropy import cosmology
cosmo = cosmology.FlatLambdaCDM(H0=70., Om0=0.3)


__all__ = ["FundamentalPlane", "Bernardi03_r"]

class FundamentalPlane(object):
    def __init__(self, a, b, c):
        """
        Define fundamental plane
        log R_0 = a log(sigma) + b log(I0) + c
        """
        self.a = a
        self.b = b
        self.c = c

    def R0(self, sigma, flux_mag, z):
        """
        Calculate effective circular size in kpc
        """
        DA = cosmo.angular_diameter_distance(z).value * 1000. / 206265.
        R0_kpc = 10**( (self.a*np.log10(sigma) \
                    - self.b/2.5*(flux_mag+2.5*np.log10(2*np.pi/DA**2)) \
                    + self.c)/(1+2*self.b))      
        return R0_kpc

    def r0(self, sigma, flux_mag, z):
        """
        Calculate effective circular size in arcsec
        """
        DA = cosmo.angular_diameter_distance(z).value * 1000. / 206265.
        return self.R0(sigma, flux_mag, z) / DA
    

Bernardi03_r = FundamentalPlane(a=1.17, b=-0.75, c=-8.022)


if __name__ == '__main__':
        
    master = Table.read('SampleZMprobaEllSub_visual.fits')

    petroflux = master['PETROFLUX'][:,4]  # nanomaggies
    petromag = 22.5 - 2.5*log10(petroflux)
    sigma = master['VDISP']
    z = master['Z_1']
    size = Bernardi03_r.R0(sigma, petromag, z)

    figure()
    scatter(log10(sigma), log10(size))
    figure()
    scatter(petromag, log10(size))
