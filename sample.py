""" Script to select sample from NSA catalog """

from pylab import *
from astropy.table import Table

rc('lines', linewidth=2)
rc('font', size=15)
rc('axes', grid=True)
rc('legend', frameon=False)

if 1:
    master = Table.read('nsa_morph.fits')
b = dict(F=0, N=1, u=2, g=3, r=4, i=5, z=6)

print 'parent      : %8i' % len(master)
# redshift
bool_Z = master['Z_1'] < 0.05
# Red sequence
colorgi = master['ABSMAG'][:,b['g']] - master['ABSMAG'][:,b['i']]
M_r = master['ABSMAG'][:,b['r']]
bool_RS = colorgi > -0.05*(M_r + 16.0) + 0.65  # from Zhu+10
# velocity dispersion
bool_VD = master['VDISP'] > 70.  # SDSS instrumental broadening
# ellipticity
ell = 1. - master['SERSIC_BA']
bool_Ell = ell < 0.6

print 'z < 0.05    : %8i' % count_nonzero(bool_Z)
print 'Red sequence: %8i' % count_nonzero(bool_RS & bool_Z)
print 'VD > 70 km/s: %8i' % count_nonzero(bool_VD & bool_RS & bool_Z)
print 'Ell < 0.6   : %8i' % count_nonzero(bool_Ell & bool_VD & bool_RS & bool_Z)
bool_OK = bool_Ell & bool_VD & bool_RS & bool_Z

# has morphology from Huertas+11
bool_MATCH = isnan(master['id']) == False
print 'has morph   : %8i' % count_nonzero(bool_MATCH & bool_OK)
bool_OKM = bool_OK & bool_MATCH

# prob dist of Sample 1
ax = figure().add_subplot(111)
bins = arange(0., 1.01, 0.1)
for k in ['probaE', 'probaEll', 'probaS0']:
    ax.hist(master[k][bool_OKM], histtype='step', bins=bins, label='%s' % k)
xlabel('prob')
ylabel('count')
legend()

bool_OKM_reject = bool_OKM & (master['probaEll'] < 0.2)
bool_OKM_select = bool_OKM & (master['probaEll'] > 0.6)
print 'Ell < 0.2   : %8i' % count_nonzero(bool_OKM_reject)
print 'Ell > 0.6   : %8i' % count_nonzero(bool_OKM_select)


# ------------
print '-'*80
print 'has morph   : %8i' % count_nonzero(bool_MATCH)

ax = figure().add_subplot(111)
bins = arange(0., 1.01, 0.1)
for k in ['probaE', 'probaEll', 'probaS0']:
    n, bins, patches = ax.hist(master[k][bool_MATCH], histtype='step', bins=bins, label='%s' % k)

bool_MATCH_reject = bool_MATCH & (master['probaEll'] < 0.2)
bool_MATCH_select = bool_MATCH & (master['probaEll'] > 0.6)

print 'Ell < 0.2   : %8i' % count_nonzero(bool_MATCH_reject)
print 'Ell > 0.6   : %8i' % count_nonzero(bool_MATCH_select)

# Find commons
print 'Ell < 0.2 common : %8i' % count_nonzero(bool_OKM_reject & bool_MATCH_reject)
print 'Ell > 0.6 common : %8i' % count_nonzero(bool_OKM_select & bool_MATCH_select)



def randTrue(boolarr, size=10):
    """ Return randomly selected index of boolarr=True """
    ind = where(boolarr)[0]
    sub = randint(len(ind), size=size)
    return ind[sub]

# names = ['IAUNAME', 'SUBIDR']

names = ['Z_1', 'RA_1', 'DEC_1']
master[names][randTrue(bool_OK, size=200)].write('rand_OK.cat', format='ascii')
master[names][randTrue(bool_MATCH_reject, size=200)].write('rand_MATCH_reject.cat', format='ascii')
master[names][randTrue(bool_MATCH_select, size=200)].write('rand_MATCH_select.cat', format='ascii')
