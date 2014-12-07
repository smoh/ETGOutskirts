# Add SIZEDEG column to the fits table for making stamp images
# Size is calculated as twice the radius at which the surface brightness is 30 mag/arcsec^2
# with a maximum of 1000pix (0.11 deg).
from pylab import *
from astropy.table import Table
from pybdfitter.sersic import bn, fn, nanomaggie2mag, sersic2d

nsa = Table.read('sample.fits')
mag_thres = 30.  # mag/arcsec^2; threshold magnitude

Re = nsa['SERSIC_TH50'] / 0.396  # pixels
ba = nsa['SERSIC_BA']
n = nsa['SERSIC_N']
Ie_avg = nsa['SERSICFLUX'][:,4] / (2.*pi*ba*(Re)**2)  # nanomaggy/pixel
Ie = Ie_avg / fn(n)  # nanomaggie/pixel
mu_e = nanomaggie2mag(Ie/.396**2)  # mag/arcsec^2

x30 = (1. + log(10)/2.5/bn(n) * (mag_thres - mu_e))**n
size_pix30 = x30*Re*2
# set maximum size to 1000pix
size_pix30[size_pix30>1000] = 1000
# pixel to degree
size_pix30 *= 0.396/3600.
size_pix30.name = 'SIZEDEG'
nsa.add_column(size_pix30)

nsa.write('sample.fits', overwrite=True)
