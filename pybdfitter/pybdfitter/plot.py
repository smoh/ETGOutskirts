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




# def plot_residual(iauname, param_fit, datadir, modeldir,
#                     savename=None, axes_title=False):
#     """
#     Plot data, model, residual
#     """

#     # load images
#     imagefits = datadir+'/deblended/%s.fits' % (iauname)
#     modelfits = modeldir+'/M%s.fits' % (iauname)
#     ivarfits = datadir+'/ivar/%s.fits' % (iauname)
#     model = fits.getdata(modelfits, ext=0)
#     data = fits.getdata(imagefits, ext=0)  # r-band
#     ivar = fits.getdata(ivarfits, ext=0)
#     residual = data - model
#     # crop images
#     if 0:
#         re = int(param_fit.Re)
#         xc, yc = int(param_fit.xc), int(param_fit.yc)
#         data = data[yc-10*re:yc+10*re, xc-10*re:xc+10*re]
#         model = model[yc-10*re:yc+10*re, xc-10*re:xc+10*re]
#         residual = residual[yc-10*re:yc+10*re, xc-10*re:xc+10*re]

#     fig = figure(figsize=(8, 7.5))
#     transData2Axes = lambda ax, p: ax.transAxes.inverted().transform(ax.transData.transform(p))

#     def add_effective_ellipse(ax):
#         from matplotlib.patches import Ellipse
#         e1 = Ellipse(xy=(param_fit.x, param_fit.y), width=2*param_fit.Re,
#                         height=2*param_fit.Re/param_fit.q, angle=rad2deg(param_fit.pa)+90.,
#                         fc='None', ec='c', lw=2)
#         e2 = Ellipse(xy=(param_fit.x, param_fit.y), width=8*param_fit.Re,
#                         height=8*param_fit.Re/param_fit.q, angle=rad2deg(param_fit.pa)+90.,
#                         fc='None', ec='c', lw=2)
#         # ax.add_artist(e)
#         ax.add_artist(e1)
#         ax.add_artist(e2)

#     # color normalization
#     norm = AsinhNorm(vmin=percentile(data.flatten(), 10), vmax=percentile(data.flatten(),90))
#     # norm = 'asinh'

#     gImages = GridSpec(1, 3)
#     gImages.update(left=0.05, right=0.95, top=0.95, bottom=0.7, wspace=0.02)
#     # data image
#     ax1 = subplot(gImages[0,0])
#     ax1, ax1_cb = showim(ax1, data, norm=norm)
#     add_effective_ellipse(ax1)

#     # best-fit model
#     ax2 = subplot(gImages[0,1], sharex=ax1, sharey=ax1)
#     ax2, ax2_cb = showim(ax2, model, norm=norm)
#     add_effective_ellipse(ax2)

#     # residual image
#     ax3 = subplot(gImages[0,2], sharex=ax1, sharey=ax1)
#     ax3, ax3_cb = showim(ax3, residual, norm=norm)
#     #showim(ax3, residual*sqrt(ivar), norm='asinh')
#     add_effective_ellipse(ax3)

#     # hide colorbar ticklabels
#     for cb in [ax1_cb, ax2_cb, ax3_cb]:
#         cb.ax.xaxis.set_ticklabels([])

#     # profiles
#     gProfile = GridSpec(4, 1)
#     gProfile.update(left=0.1, right=0.9, top=0.65, bottom=0.1)
#     ax4 = subplot(gProfile[:3,0])
#     # divider = make_axes_locatable(ax4)
#     # ax_res = divider.append_axes("bottom", size="30%", pad=0.1, sharex=ax4)
#     ax_res = subplot(gProfile[3,0], sharex=ax4)
#     fig.add_axes(ax_res)

#     pData = getProfile(
#             data, param_fit.Re, param_fit.q, param_fit.pa, param_fit.x, param_fit.y)
#     pModel = getProfile(
#             model, param_fit.Re, param_fit.q, param_fit.pa, param_fit.x, param_fit.y)
#     pRes = getProfile(
#             residual, param_fit.Re, param_fit.q, param_fit.pa, param_fit.x, param_fit.y)
#     ivar[isnan(ivar)] = 1.
#     pVar = getProfile(
#             1./ivar, param_fit.Re, param_fit.q, param_fit.pa, param_fit.x, param_fit.y)
#     pRedChi = getProfile(
#             residual**2*ivar,
#             param_fit.Re, param_fit.q, param_fit.pa, param_fit.x, param_fit.y,
#             limit=-1)

#     rad = pData['rad'] / param_fit.Re
#     sbData = -2.5*log10(pData['sb']/0.396**2) + 22.5
#     sbModel = -2.5*log10(pModel['sb']/0.396**2) + 22.5
#     sbSigma = np.abs(2.5*sqrt(pVar['totalflux']) / pData['totalflux'] / log(10.) * sqrt(pVar['area']*0.396**2))

#     ax4.plot(rad, sbData, 'o-', color='0.5')
#     # ax4.errorbar(rad, sbData, yerr=sbSigma, fmt='ko-')
#     good = sbSigma < 2  # to manage dynamical range
#     ax4.fill_between(rad, sbData-sbSigma, sbData+sbSigma, where=good,
#                      edgecolor='None', facecolor='0.85')
#     ax4.plot(rad, sbModel, 'r-', lw=2)
#     # ax4.axvline(rad[good][-1], ls='--', c='k')  # radius at which sigma(mag/arsec^2) reaches 2
#     # ax_res.errorbar(rad, sbData-sbModel, yerr=sbSigma, fmt='bo-')
#     ax_res.errorbar(rad, (sbData-sbModel)/sbSigma, fmt='bo-')
#     ax_res.axhline(0, ls=':', c='k')
#     # ax_res.plot(rad, pRes['mnflux'], 'ko-')
#     # ax_res.plot(rad, pData['mnflux']-pModel['mnflux'], 'r-')
#     mpltools.hide_tick_labels(ax4, 'x')
#     ax4.invert_yaxis()
#     ax4.minorticks_on()
#     ax4.set_ylabel(r"$\mu$ (mag/arcsec$^2$)")
#     ax_res.minorticks_on()
#     ax_res.set_xlabel("$R/R_e$")
#     ax_res.set_ylim(-5, 5)

#     ax_chi = twinx(ax4)  #.twinx()
#     radchi = pRedChi['rad']/param_fit.Re
#     ax_chi.plot(radchi, pRedChi['mnflux'], 's-', c='b', mec='b', mfc='None')
#     # cumulative reduced chi-squre
#     cumChi = cumsum(pRedChi['totalflux'])
#     cumN = cumsum(pRedChi['area'])
#     cumRedChi = cumChi / cumN
#     ax_chi.plot(radchi, cumRedChi, '.-', c='b')
#     ax_chi.set_yscale('log')
#     ax_chi.set_xlim(0, 8)
#     ax_chi.axhline(1., ls='--', c='b')

    
#     subplots_adjust(left=0.03, right=0.95, top=0.95, bottom=0.05, wspace=0.01)
#     mpltools.hide_tick_labels(ax2, 'y')
#     mpltools.hide_tick_labels(ax3, 'y')

#     if 0:
#         for ax in [ax1, ax2, ax3]:
#             ax.set_xlim(param_fit.xc - param_fit.Re*5, param_fit.xc + param_fit.Re*5)
#             ax.set_ylim(param_fit.yc - param_fit.Re*5, param_fit.yc + param_fit.Re*5)

#     if savename:
#         fig.savefig(savename)
#         close()

#     return vars()


# def test_plot_residual():
#     if 0:
#         result = Table.read('sdss_psf/ser/deblended/RAWFIT00000.00048.fits')
#         datadir = 'sdss_psf/ser/data/'
#         modeldir = 'sdss_psf/ser/deblended/models'
#         gal = result[9]
#         v = plot_residual(gal['NAME'], Sersic(gal['FIT_SER']), datadir, modeldir)
#     if 0:
#         result = Table.read('sdss_psf/ser/deblended/RAWFIT00000.00048.fits')
#         datadir = 'sdss_psf/ser/data/'
#         modeldir = 'sdss_psf/ser/deblended/models'
#         for i, gal in enumerate(result):
#             print i
#             v = plot_residual(gal['NAME'], Sersic(gal['FIT_SER']), datadir, modeldir,
#                             savename='sdss_psf/ser/deblended/plotProfile/%s.png' % gal['NAME'])
#     if 1:
#         result = Table.read('sdss_psf/dvc/deblended/RAWFIT00000.00048.fits')
#         datadir = 'sdss_psf/dvc/data/'
#         modeldir = 'sdss_psf/dvc/deblended/models'
#         for i, gal in enumerate(result):
#             print i
#             v = plot_residual(gal['NAME'], Sersic(gal['FIT_DVC']), datadir, modeldir,
#                             savename='sdss_psf/dvc/deblended/plotProfile/%s.png' % gal['NAME'])
#     return v


# def main():
#     """
#     Save image/model/residual plots
#     """
#     import argparse
#     parser = argparse.ArgumentParser(
#                 description="""Save image/model/residual plots""")
#     parser.add_argument("input", type=str, help="part of NSA catalog")
#     parser.add_argument("result", type=str, help="fit result table RAWXXX")
#     parser.add_argument("profile", type=str, help="name of profile")
#     parser.add_argument("datadir", type=str, help="input image directory")
#     parser.add_argument("modeldir", type=str, help="model image directory")
#     parser.add_argument("outdir", type=str, help="output directory")
#     parser.add_argument("-r", "--range", nargs=2, type=int,
#                         help="zero-based range of row index, e.g., -r 0 5")
#     args = parser.parse_args()

#     mpltools.mkdir_p(args.outdir)

#     inputt = Table.read(args.input)
#     result = Table.read(args.result)
#     datadir, modeldir, profile, outdir = \
#          args.datadir, args.modeldir, args.profile, args.outdir
#     row_start = args.range[0] if args.range else 0
#     row_end = args.range[1] + 1 if args.range else len(result)

#     for index in range(row_start, row_end):
#         gal = result[index]
#         param_fit = Sersic(gal['FIT_%s' % (profile)])
#         plot_residual(gal['NAME'], param_fit, datadir, modeldir,
#             savename=outdir+'/%s.png'%(gal['NAME']))


# if __name__ == '__main__':
#     main()