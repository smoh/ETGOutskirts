from pylab import *
from astropy.table import Table

import mpltools
mpltools.set_pub_single()

master = Table.read('../SampleZMprobaEllSub_visual.fits')
dvc = Table.read('dvc/deblended/RAWFIT00000.00048.fits')
ser = Table.read('ser/deblended/RAWFIT00000.00048.fits')

def chisq():
    # scatter(dvc['CHISQ_DVC'], ser['CHISQ_SER']/dvc['CHISQ_DVC'])

    def scatter_hist(ax, x, y):

        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        ax_sub = divider.append_axes('right', size='30%', pad=.2, sharey=ax)
        mpltools.hide_tick_labels(ax_sub, 'y')

        ax.scatter(x, y)

        ax_sub.hist(y, orientation='horizontal')

        return ax, ax_sub

    ax, ax_sub = scatter_hist(gca(), dvc['CHISQ_DVC'], ser['CHISQ_SER']/dvc['CHISQ_DVC'])

    ax.axhline(1., ls='--', c='k')
    ax.set_xlabel('CHISQ_DVC')
    ax.set_ylabel('CHISQ_SER / CHISQ_DVC')

    ax_sub.set_xlim(0, 30)


def plot_identity(ax, **kwargs):

    x1, x2 = ax.get_xlim()
    y1, y2 = ax.get_ylim()
    x = max(x1, y1)
    y = min(x2, y2)
    v = linspace(x, y, 100)
    ax.plot(v, v, **kwargs)


def Reff():
    pser = ser['FIT_SER'].data
    pdvc = dvc['FIT_DVC']
    
    scatter(pser[:,1], pdvc[:,1])
    plot(linspace(0, 50, 100), linspace(0, 50, 100))

    xlim(0, 120)
    ylim(0, 50)
    xlabel('Reff_SER (pix)')
    ylabel('Reff_DVC (pix)')


def FP_Bernardi03(sigma, I0):
    return 10**(1.17*log10(sigma) - 0.75*log10(I0) - 8.022)

def bn_Sersic(n):
    b = n*exp(0.6950-0.1789/n)
    return b

def totalflux_Sersic(n, Ie, Reff, q_ba):
    """ total flux of one Sersic model """
    from scipy.special import gamma
    b = bn_Sersic(n)
    flux = 2*pi*n * Reff**2 * Ie * exp(b) / b**(2.*n) * gamma(2.*n)
    return flux


# def fp():
if 1:
    from astropy import cosmology
    cosmo = cosmology.FlatLambdaCDM(H0=70., Om0=0.3)

    z = master['Z_1'].data
    n = dvc['FIT_DVC'].data[:,2]
    Ie = dvc['FIT_DVC'].data[:,0]
    reff_pix = dvc['FIT_DVC'].data[:,1]
    q_ba = dvc['FIT_DVC'].data[:,3]
    
    DA = cosmo.angular_diameter_distance(z).value * 1000. / 206265.  # kpc/arcsec
    R0 = dvc['FIT_DVC'].data[:,1] * 0.396 * DA * q_ba  # kpc

    sigma = master['VDISP'].data

    flux = totalflux_Sersic(n, Ie, reff_pix, q_ba)  # nanomaggies
    area = pi*(reff_pix * q_ba)**2
    I0 = flux*1e-9 / (area * 0.396**2)  # Jy / arcsec^2
    lhs = 1.17 * log10(sigma) - 0.75 * log10(I0) - 8.022
    mu0 = -2.5*log10(I0)
    lhs1 = log10(sigma) + 0.20 * (mu0-20.09)
    scatter(lhs1, log10(R0))

    # scatter(log10(sigma), log10(R0))

    # x = arange(-10, -8.5, 0.1)
    # plot(x-8.022, x-8.022, 'k--')
    
    xlabel('a log(sigma) + b log(I0)')
    ylabel('log Re')

