""" Plotting fitting results """

import pylab as pl
import numpy as np
from astropy.table import Table
from astropy.io import fits
from asinh_norm import AsinhNorm
from matplotlib.colors import LogNorm
import math as m
from mpl_toolkits.axes_grid1 import make_axes_locatable

from sersic import Sersic, NSersic, sersic2d, nanomaggie2mag
import mpltools
from itertools import cycle

pl.seterr(all='warn')
pl.rc('lines', lw=2)


def add_ellipse(ax, x, y, rx, ry, pa, **kwargs):
    """
    Add an ellipse to the axes

    ax : matplotlib axes instance
    x, y : center
    rx, ry : semi- major, minor axis
    pa : radian from x axis counter-clockwise
    kwargs : passed on to matplotlib.patches.Ellipse
    """
    # default
    if not 'fc' in kwargs:
        kwargs['fc'] = 'None'
    if not 'ec' in kwargs:
        kwargs['ec'] = 'c'
    if not 'lw' in kwargs:
        kwargs['lw'] = 1.

    from matplotlib.patches import Ellipse
    e = Ellipse(xy=(x, y), width=2*rx, height=2*ry, angle=pl.rad2deg(pa),
                **kwargs)
    ax.add_artist(e)


class Ellipse:
    def __init__(self, r, q, x0, y0, phi):
        """
        Define an ellipse

        r : semi-major axis
        q : minor/major axis ratio
        x0 : center x coordinate
        y0 : center y coordinate
        phi : major axis position angle.
               counter-clockwise from x-axis in radian
        """
        self.r = r
        self.x0 = x0
        self.y0 = y0
        self.phi = phi
        self.q = q

    def recenter(self, x, y):
        """
        Return (x, y) coordinate relative to the center of the ellipse with
        x-axis along the major axis
        """
        return ((x-self.x0)*np.cos(self.phi) + (y-self.y0)*np.sin(self.phi),
                (y-self.y0)*np.cos(self.phi) - (x-self.x0)*np.sin(self.phi))

    def is_inside(self, x, y):
        """
        Retrun true if x,y point is in ellipse
        """
        (xt, yt) = self.recenter(x, y)
        inside = (xt < self.r) & (yt < self.r)
        small = np.where(inside)

        rads = np.sqrt(xt[small] * xt[small] +
                       yt[small] * yt[small] / (self.q * self.q))
        inside[small] = rads < self.r
        return inside


def calc_annulus(image, r_in, r_out, q, phi, x0, y0, mask=None):
    """
    Determine the sum, mean, variance of values within a elliptical annulus

    image : 2-d numpy array
    r_in : inner semi-major axis
    r_out : outer semi-major axis
    q : major/minor axis ratio
    phi : position angle, radian from x axis counter-clockwise
    x0, y0 : center of ellipse
    mask : mask array with the same shape as image
            True for masked

    Returns (total, mean, variance, area)
    area is -1 if total = 0
    """
    Ny, Nx = pl.shape(image)
    if mask is None:
        mask = pl.array([[False]*Nx]*Ny)
    x, y = pl.meshgrid(pl.arange(Nx), pl.arange(Ny))
    ell1 = ~(Ellipse(r_in, q, x0, y0, phi).is_inside(x, y))
    ell2 = Ellipse(r_out, q, x0, y0, phi).is_inside(x, y)
    annulus = ell1 & ell2 & ~mask
    if not annulus.any():
        return (0.0, 0.0, 0.0, -1.0)
    totalflux, meanflux, varflux, area = 0., 0., 0., 0.
    totalflux = sum(image[annulus])
    area = image[annulus].size
    meanflux = pl.mean(image[annulus])
    varflux = pl.var(image[annulus])
    if totalflux == 0.0:
        area = -1.0
    return (totalflux, meanflux, varflux, area)


def get_profile(image, q, phi, x0, y0, step=5, limit=100, mask=None):
    """
    Get elliptical annulus profile
    """
    r = 0.
    numstep = int(limit/step)
    mult_factor = 1.
    profile = np.recarray((numstep,),
                          dtype=[('radius', float),
                                 ('totalflux', float),
                                 ('mnflux', float),
                                 ('stdflux', float),
                                 ('sb', float),
                                 ('area', float)])
    for i in range(numstep):
        profile[i].radius = r + 0.5*step
        tf, mf, vf, area = calc_annulus(image, r, r+step, q, phi, x0, y0,
                                        mask=mask)
        profile[i].mnflux = mf
        profile[i].stdflux = m.sqrt(vf)
        profile[i].sb = tf/area
        profile[i].area = area
        profile[i].totalflux = tf
        r = r + step
        step *= mult_factor

    return profile


def showim(ax, image, norm=None, **kwargs):
    """
    Show single image with colorbar on the bottom

    ax : matplotlib axes instance
    image : 2-d image array

    Keywords
    --------
    norm : color normalization. Either Normalize instance or one of following
        'asinh' : arcsinh
        'log' : LogNorm
        none : default

    Returns
    tuple of image axes, colorbar axes
    """

    divider = make_axes_locatable(ax)
    ax_cb = divider.append_axes("top", size="5%", pad=0.05)
    cf = ax.get_figure()
    cf.add_axes(ax_cb)

    if not isinstance(norm, pl.Normalize):
        if norm is 'asinh':
            norm = AsinhNorm(
                # vmin=image[~isnan(image)].min(),
                vmin=pl.percentile(image[~pl.isnan(image)].flatten(), 5),
                vmax=pl.percentile(image[~pl.isnan(image)].flatten(), 90)*20)
        if norm is 'log':
            norm = LogNorm()
        if norm is None:
            pass
    im = ax.imshow(image, cmap=pl.cm.jet,
                   norm=norm,
                   interpolation='nearest',
                   aspect='auto', origin='lower', **kwargs)
    # ax.xaxis.tick_top()
    try:
        cb = pl.colorbar(im, cax=ax_cb, orientation='horizontal')
        ax_cb.xaxis.tick_top()
        return ax, cb
    except:
        pass


class DataContainer(object):
    """ DataContainer class """
    def __init__(self, img, ivar, model, p):
        """
        Load data, model and residual images
        """
        self.img = fits.getdata(img)
        self.ivar = fits.getdata(ivar)
        self.model = fits.getdata(model)
        self.residual = self.img - self.model
        self.param = NSersic(p)


class Plotter(object):
    """plotting fitting result"""
    def __init__(self, output, datadir='data', modeldir='models'):
        """
        Initialize Plotter class

        output : str, output filename (RAWFITXXX.fits)
        datadir : data directory
        modeldir : models directory
        """
        if type(output) == str:
            self.result = Table.read(output)
        else:
            self.result = output
        self.datadir = datadir
        self.modeldir = modeldir
        self._get_profile_name()

    def _get_profile_name(self):
        # find out profile name
        fitcol = [k for i, k in enumerate(self.result.colnames) if 'FIT_' in k]
        assert len(fitcol) == 1, "More than one profile found!"
        self.profile_name = fitcol[0].split('FIT_')[1]

    def __getitem__(self, index):
        return DataContainer(
            self.datadir+'/deblended/%s.fits' % (self.result[index]['NAME']),
            self.datadir+'/ivar/%s.fits' % (self.result[index]['NAME']),
            self.modeldir+'/M%s.fits' % (self.result[index]['NAME']),
            self.result[index]['FIT_'+self.profile_name])

    def show_images(self, index):
        """
        Show data, model, and residual image
        """
        img = self[index].img
        model = self[index].model
        residual = img - model
        sersic = NSersic(self.result[index]['FIT_'+self.profile_name])

        fig = pl.figure(figsize=(12, 4))
        gImages = pl.GridSpec(1, 3)
        gImages.update(
            left=0.05, right=0.97, top=0.9, bottom=0.1, wspace=0.02)
        # color normalization
        norm = AsinhNorm(
            vmin=pl.percentile(img.flatten(), 10),
            vmax=pl.percentile(img.flatten(), 90))
        # data image
        ax1 = pl.subplot(gImages[0, 0])
        ax1, ax1_cb = showim(ax1, img, norm=norm)
        # best-fit model
        ax2 = pl.subplot(gImages[0, 1], sharex=ax1, sharey=ax1)
        ax2, ax2_cb = showim(ax2, model, norm=norm)
        # residual image
        ax3 = pl.subplot(gImages[0, 2], sharex=ax1, sharey=ax1)
        ax3, ax3_cb = showim(ax3, residual, norm=norm)
        # highlight effective radius of each component
        colors = cycle(['#3CF5C1', '#F53C70'])
        for i in range(sersic.nprofiles):
            add_ellipse(ax3, sersic[i].x, sersic[i].y, sersic[i].Re,
                        sersic[i].Re*sersic[i].q, sersic[i].pa,
                        ec=colors.next())
        # hide ticklabels
        mpltools.hide_tick_labels(ax2, 'y')
        mpltools.hide_tick_labels(ax3, 'y')
        return fig

    def show_profile(self, index):
        """
        Show elliptical annulus profile
        """
        sersic = NSersic(self.result[index]['FIT_'+self.profile_name])

        fig = pl.figure()
        gProfile = pl.GridSpec(4, 1)
        ax = pl.subplot(gProfile[:3, 0])
        ax_res = pl.subplot(gProfile[3, 0], sharex=ax)
        fig.add_axes(ax_res)

        #TODO: manage compute time with step (divide limit by some nstep?)
        limit = 2.*max([sersic[i].Re for i in range(sersic.nprofiles)])
        step = int(limit/20.)

        # deal with <=zero variance
        ivar = self[index].ivar.copy()
        mask = ivar <= 1e-30
        ivar[mask] = 1.
        v = 1./ivar

        azimuth = lambda x: get_profile(x, 1., 0., sersic[0].x, sersic[0].y,
                                        step=step, limit=limit, mask=mask)

        pVar = azimuth(v)
        # sigma = pl.sqrt(pVar.totalflux)/pVar.area  #pl.sqrt(pVar.area)
        sigma = pl.sqrt(pVar.mnflux)
        pImg = azimuth(self[index].img)
        # sigma = pImg.stdflux
        pModel = azimuth(self[index].model)
        good = (sigma < pImg.mnflux) & (pImg.mnflux > 0)

        if 1:
            ax.plot(pImg.radius, pImg.mnflux, c='k')
            # ax.fill_between(pImg.radius[good],
            #                  (pImg.mnflux-sigma)[good],
            #                  (pImg.mnflux+sigma)[good],
            #                  edgecolor='None', facecolor='0.5')
            ax.plot(pModel.radius, pModel.mnflux, 'r-')
            ax.set_yscale('log')
            ax.set_ylabel('nnmg/pixel')
        if 0:
            radius = pImg.radius[good]
            sbImg = nanomaggie2mag(pImg.mnflux[good]/0.396**2)
            sbModel = nanomaggie2mag(pModel.mnflux[good]/0.396**2)
            sbSigma = 2.5*sigma/pl.log(10.)/pImg.mnflux[good]
            ax.plot(pImg.radius, sbImg, c='k')
            ax.fill_between(radius,
                            sbImg-sbSigma,
                            sbImg+sbSigma,
                            edgecolor='None', facecolor='0.5')
            ax.plot(pModel.radius, sbModel, 'r-')
            ax.set_ylim(max(sbImg), min(sbImg))
            ax.set_ylabel('mag/arcsec$^2$')
        # for multi component fit
        # these subcomponents are NOT psf-convolved!
        if 1:
            if sersic.nprofiles > 1:
                x, y = pl.meshgrid(pl.arange(self[index].img.shape[1]),
                                   pl.arange(self[index].img.shape[0]))
                complist = []
                for i in range(sersic.nprofiles):
                    complist.append(sersic2d(x, y, sersic[i].p))
                clist   = cycle(['g', 'm'])
                for im in complist:
                    p = azimuth(im)
                    ax.plot(p.radius, p.mnflux,
                            c=clist.next())
        # vertical lines indicating N X Re
        for i in range(sersic.nprofiles):
            for n in pl.arange(1, 1.1):
                ax.axvline(sersic[i].Re*n, ls='--')
        ax.set_xlim([0, pImg.radius[-1]])

        # residual
        ax_res.plot(pImg.radius, (pImg.mnflux - pModel.mnflux)/sigma)

        return vars()




