from pylab import *
from astropy.table import Table

from mpltools import get_handle as gh

def plotCoeffsLiterature():
    """
    Plot Fundamental Plane coefficients a, b, c from the literature
    """
    t = Table.read('fp_coeff.csv', delimiter=',', format='ascii')

    def str2float(s):
        s1 = s.replace('\xe2\x88\x92','-')
        try:
            out = float(s1)
        except ValueError:
            out = 0
        return out

    a = array([str2float(x.split(' ')[0]) for x in t['a']])
    a_err = array([str2float(x.split(' ')[-1]) for x in t['a']])
    b = array([str2float(x.split(' ')[0]) for x in t['b']])
    b_err = array([str2float(x.split(' ')[-1]) for x in t['b']])
    c = array([str2float(x.split(' ')[0]) for x in t['c']])
    c_err = array([str2float(x.split(' ')[-1]) for x in t['c']])

    # hard fix
    c_err[t['Authors']=='Hyde & Bernardi (2009)'] = 0.

    i_band = array([ ['u','g','r','i','z'].index(x) for x in t['Band'] ])
    cdict = {'Bernardi et al. (2003c)':'k', 'Hyde & Bernardi (2009)':'b', 'Saulder et al. (2013)':'r',
             'La Barbera et al. (2010a)':'g', 'La Barbera et al. (2008)':'g', 'Jorgensen et al. (1996)':'m'}
    lsdict = {'direct ML':'-', 'direct R':'-', 'orthogonal ML':':', 'orthogonal R':':'}

    fig = figure(figsize=(8,6))
    ax1 = fig.add_subplot(221)
    ax1.set_title('a')
    ax2 = fig.add_subplot(222)
    ax2.set_title('b')
    ax3 = fig.add_subplot(223)
    ax3.set_title('c')
    ax3.set_ylim(-9.5, -7.0)
    for ind, src in enumerate(sorted(cdict.keys())):
        for fit in lsdict.keys():
            sub = (t['Authors']==src) & (t['Type of fit'] == fit)
            x = i_band[sub]

            for coeff, coeff_err, ax in zip([a,b,c],[a_err,b_err,c_err],[ax1,ax2,ax3]):
            
                val = coeff[sub][x.argsort()]
                err = coeff_err[sub][x.argsort()]

                ax.errorbar(x + (ind-3)*0.1, val, yerr=err, fmt='o-', c=cdict[src], ls=lsdict[fit])
        
    for ax in fig.get_axes():
        ax.set_xticks([0,1,2,3,4])
        ax.set_xticklabels(['u','g','r','i','z'])
        ax.set_xlim(-0.5, 4.5)

    # legend
    handles = [gh(c='k', ls=':'), gh(c='k', ls='-')]
    labels = ['orthogonal' ,'direct']
    for src in cdict.keys():
        handles.append(gh(c=cdict[src], ls='-'))
        labels.append(src)
    leg = figlegend(handles, labels, 'lower right')
    [l.set_linewidth(2) for l in leg.get_lines()]


    figtext(0.5, 0.93, 'Fundamental Plane', ha='center', va='center')


from mpl_toolkits.mplot3d import Axes3D

fig = figure()
ax = fig.add_subplot(111, projection='3d')

# FP coefficients
# log10(R0) = a * log10(sigma0) + b * log10(I0) + c

X = arange(2., 2.6, 0.01)  # log sig0
Y = arange(-7, -10, -0.01)  # log I0

XX, YY = meshgrid(X, Y)
Z = 1.5*XX - 0.8*YY -8.7

ax.plot_surface(XX, YY, Z)