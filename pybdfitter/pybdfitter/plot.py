""" Script to plot fitting results """

from pylab import *
from astropy.table import Table
from astropy.io import fits
import aplpy
from asinh_norm import AsinhNorm
from matplotlib.colors import LogNorm
import numpy as np
import math as m
from mpl_toolkits.axes_grid1 import make_axes_locatable

from sersic import Sersic
import mpltools

rc('savefig', dpi=80)


def add_ellipse(ax, x, y, rx, ry, pa,**kwargs):
    """
    Add an ellipse to the axes

    ax : axes
    x, y : center
    rx, ry : semi- major, minor axis
    pa : radian from x axis counter-clockwise
    kwargs : passed on to matplotlib.patches.Ellipse
    """
    # default
    if not kwargs.has_key('fc'):
        kwargs['fc'] = 'None'
    if not kwargs.has_key('ec'):
        kwargs['ec'] = 'c'
    if not kwargs.has_key('lw'):
        kwargs['lw'] = 2

    from matplotlib.patches import Ellipse
    e = Ellipse(xy=(x, y), width=2*rx, height=2*ry, angle=rad2deg(pa)+90.,
                **kwargs)
    ax.add_artist(e)


def azimuthalAverage(image, center=None, stddev=False, returnradii=False, return_nr=False, 
        binsize=0.5, weights=None, steps=False, interpnan=False, left=None, right=None,
        mask=None ):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fractional pixels).
    stddev - if specified, return the azimuthal standard deviation instead of the average
    returnradii - if specified, return (radii_array,radial_profile)
    return_nr   - if specified, return number of pixels per radius *and* radius
    binsize - size of the averaging bin.  Can lead to strange results if
        non-binsize factors are used to specify the center and the binsize is
        too large
    weights - can do a weighted average instead of a simple average if this keyword parameter
        is set.  weights.shape must = image.shape.  weighted stddev is undefined, so don't
        set weights and stddev.
    steps - if specified, will return a double-length bin array and radial
        profile so you can plot a step-form radial profile (which more accurately
        represents what's going on)
    interpnan - Interpolate over NAN values, i.e. bins where there is no data?
        left,right - passed to interpnan; they set the extrapolated values
    mask - can supply a mask (boolean array same size as image with True for OK and False for not)
        to average over only select data.

    If a bin contains NO DATA, it will have a NAN value because of the
    divide-by-sum-of-weights component.  I think this is a useful way to denote
    lack of data, but users let me know if an alternative is prefered...
    
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if center is None:
        center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    if weights is None:
        weights = np.ones(image.shape)
    elif stddev:
        raise ValueError("Weighted standard deviation is not defined.")

    if mask is None:
        mask = np.ones(image.shape,dtype='bool')
    # obsolete elif len(mask.shape) > 1:
    # obsolete     mask = mask.ravel()

    # the 'bins' as initially defined are lower/upper bounds for each bin
    # so that values will be in [lower,upper)  
    nbins = int(np.round(r.max() / binsize)+1)
    maxbin = nbins * binsize
    bins = np.linspace(0,maxbin,nbins+1)
    # but we're probably more interested in the bin centers than their left or right sides...
    bin_centers = (bins[1:]+bins[:-1])/2.0

    # how many per bin (i.e., histogram)?
    # there are never any in bin 0, because the lowest index returned by digitize is 1
    #nr = np.bincount(whichbin)[1:]
    nr = np.histogram(r,bins)[0]

    # recall that bins are from 1 to nbins (which is expressed in array terms by arange(nbins)+1 or xrange(1,nbins+1) )
    # radial_prof.shape = bin_centers.shape
    if stddev:
        # Find out which radial bin each point in the map belongs to
        whichbin = np.digitize(r.flat,bins)
        # This method is still very slow; is there a trick to do this with histograms? 
        radial_prof = np.array([image.flat[mask.flat*(whichbin==b)].std() for b in xrange(1,nbins+1)])
    else: 
        radial_prof = np.histogram(r, bins, weights=(image*weights*mask))[0] / np.histogram(r, bins, weights=(mask*weights))[0]

    if interpnan:
        radial_prof = np.interp(bin_centers,bin_centers[radial_prof==radial_prof],radial_prof[radial_prof==radial_prof],left=left,right=right)

    if steps:
        xarr = np.array(zip(bins[:-1],bins[1:])).ravel() 
        yarr = np.array(zip(radial_prof,radial_prof)).ravel() 
        return xarr,yarr
    elif returnradii: 
        return bin_centers,radial_prof
    elif return_nr:
        return nr,bin_centers,radial_prof
    else:
        return radial_prof


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

    def is_inside( self, x, y ):
        """
        Retrun true if x,y point is in ellipse
        """
        (xt,yt)=self.recenter(x,y)
        inside=(xt<self.r)&(yt<self.r)
        small=np.where(inside)

        rads=np.sqrt(xt[small]*xt[small] + yt[small]*yt[small]/(self.q*self.q))
        inside[small] = rads < self.r
        return inside


def FluxAnnulus(image, r_in, r_out, q, phi, x0, y0 ):
    """
    determine the total flux within an annulus
    also return the mean and the variance
    """
    #repixel by 1
    rp= 1.
    #image_64 = repixel(np.asarray(image), 1/rp )
    image_64 = image
    norm = np.sum(image)/np.sum(image_64)

    #image_64 = image
    x,y = np.meshgrid( np.arange(np.shape(image_64)[1]*1.0),
               np.arange(np.shape(image_64)[0]*1.0) )
    ell1 = ~(Ellipse(r_in*rp, q, rp*x0, rp*y0, phi ).is_inside(x,y))
    ell2 = Ellipse(r_out*rp, q, rp*x0, rp*y0, phi).is_inside(x,y)
    annulus = ell1 & ell2
    if ~annulus.any():
        return (0.0,0.0,0.0,-1.0)
    totalFlux,meanFlux,varFlux,area=(0.0,0.0,0.,0.)
    totalFlux = np.sum(image_64[annulus])*norm
    area = (image_64[annulus].size/(rp*rp))
    meanFlux = np.mean(image_64[annulus]*norm)
    varFlux = np.var(image_64[annulus]*norm)
    if(totalFlux == 0.0): area = -1.0
    return (totalFlux, meanFlux, varFlux, area)
               

def getProfile(image, reff, q, phi, x0, y0, limit=5):
    """
    
    limit : the limit of radial profile in unit of reff
    """

    # step = min(max(0.05*reff,0.8),1.5)
    # r=0.0
    # numstep=30
    # high = 4.*reff
    # mult_factor = 10.**(np.log10(high-step)/numstep)

    step = 5
    r = 0.0
    # numstep = int(floor( min(np.abs([x0, image.shape[1]-reff, y0, image.shape[0]-reff]))/step))
    if limit > 0:
        numstep = int(limit*reff/step)
    elif limit == -1:
        rmax = (min(image.shape)-50)/2.
        numstep = int(rmax/step)
    mult_factor = 1.
 
    profile = np.recarray((numstep,),dtype=[('rad',float),
                        ('totalflux', float),
                        ('mnflux',float),
                        ('stdflux',float),
                        ('sb', float),
                        ('area', float)])
    for i in range(numstep):
        profile[i].rad = r + 0.5*step
        (tf, mf, vf, area) = FluxAnnulus(image, r, r+step, q, phi, x0, y0)
        profile[i].mnflux = mf
        profile[i].stdflux = m.sqrt(vf)
        profile[i].sb = tf/area
        profile[i].area = area
        profile[i].totalflux = tf
        r = r + step
        step *= mult_factor

    return profile
    


def showim(ax, image, norm=norm):
    """ show single image """

    divider = make_axes_locatable(ax)
    ax_cb = divider.append_axes("bottom", size="5%", pad=0.05)
    cf = ax.get_figure()
    cf.add_axes(ax_cb)

    if not isinstance(norm, Normalize):
        if norm is 'asinh':
            norm = AsinhNorm(
                    # vmin=image[~isnan(image)].min(),
                    vmin=percentile(image[~isnan(image)].flatten(), 5),
                    vmax=percentile(image[~isnan(image)].flatten(), 90)*20)
        if norm is 'log':
            norm = LogNorm()
        if norm is None:
            pass

    im = ax.imshow(image, cmap=cm.gray_r, 
            norm=norm,
            interpolation='nearest',
            aspect='auto')
    ax.xaxis.tick_top()
    try:
        cb = colorbar(im, cax=ax_cb, orientation='horizontal')
        ax_cb.xaxis.tick_bottom()
        return ax, cb
    except:
        pass


def plot_residual(iauname, param_fit, datadir, modeldir,
                    savename=None, axes_title=False):
    """
    Plot data, model, residual
    """

    # load images
    imagefits = datadir+'/deblended/%s.fits' % (iauname)
    modelfits = modeldir+'/M%s.fits' % (iauname)
    ivarfits = datadir+'/ivar/%s.fits' % (iauname)
    model = fits.getdata(modelfits, ext=0)
    data = fits.getdata(imagefits, ext=0)  # r-band
    ivar = fits.getdata(ivarfits, ext=0)
    residual = data - model
    # crop images
    if 0:
        re = int(param_fit.Re)
        xc, yc = int(param_fit.xc), int(param_fit.yc)
        data = data[yc-10*re:yc+10*re, xc-10*re:xc+10*re]
        model = model[yc-10*re:yc+10*re, xc-10*re:xc+10*re]
        residual = residual[yc-10*re:yc+10*re, xc-10*re:xc+10*re]

    fig = figure(figsize=(8, 7.5))
    transData2Axes = lambda ax, p: ax.transAxes.inverted().transform(ax.transData.transform(p))

    def add_effective_ellipse(ax):
        from matplotlib.patches import Ellipse
        e1 = Ellipse(xy=(param_fit.x, param_fit.y), width=2*param_fit.Re,
                        height=2*param_fit.Re/param_fit.q, angle=rad2deg(param_fit.pa)+90.,
                        fc='None', ec='c', lw=2)
        e2 = Ellipse(xy=(param_fit.x, param_fit.y), width=8*param_fit.Re,
                        height=8*param_fit.Re/param_fit.q, angle=rad2deg(param_fit.pa)+90.,
                        fc='None', ec='c', lw=2)
        # ax.add_artist(e)
        ax.add_artist(e1)
        ax.add_artist(e2)

    # color normalization
    norm = AsinhNorm(vmin=percentile(data.flatten(), 10), vmax=percentile(data.flatten(),90))
    # norm = 'asinh'

    gImages = GridSpec(1, 3)
    gImages.update(left=0.05, right=0.95, top=0.95, bottom=0.7, wspace=0.02)
    # data image
    ax1 = subplot(gImages[0,0])
    ax1, ax1_cb = showim(ax1, data, norm=norm)
    add_effective_ellipse(ax1)

    # best-fit model
    ax2 = subplot(gImages[0,1], sharex=ax1, sharey=ax1)
    ax2, ax2_cb = showim(ax2, model, norm=norm)
    add_effective_ellipse(ax2)

    # residual image
    ax3 = subplot(gImages[0,2], sharex=ax1, sharey=ax1)
    ax3, ax3_cb = showim(ax3, residual, norm=norm)
    #showim(ax3, residual*sqrt(ivar), norm='asinh')
    add_effective_ellipse(ax3)

    # hide colorbar ticklabels
    for cb in [ax1_cb, ax2_cb, ax3_cb]:
        cb.ax.xaxis.set_ticklabels([])

    # profiles
    gProfile = GridSpec(4, 1)
    gProfile.update(left=0.1, right=0.9, top=0.65, bottom=0.1)
    ax4 = subplot(gProfile[:3,0])
    # divider = make_axes_locatable(ax4)
    # ax_res = divider.append_axes("bottom", size="30%", pad=0.1, sharex=ax4)
    ax_res = subplot(gProfile[3,0], sharex=ax4)
    fig.add_axes(ax_res)

    pData = getProfile(
            data, param_fit.Re, param_fit.q, param_fit.pa, param_fit.x, param_fit.y)
    pModel = getProfile(
            model, param_fit.Re, param_fit.q, param_fit.pa, param_fit.x, param_fit.y)
    pRes = getProfile(
            residual, param_fit.Re, param_fit.q, param_fit.pa, param_fit.x, param_fit.y)
    ivar[isnan(ivar)] = 1.
    pVar = getProfile(
            1./ivar, param_fit.Re, param_fit.q, param_fit.pa, param_fit.x, param_fit.y)
    pRedChi = getProfile(
            residual**2*ivar,
            param_fit.Re, param_fit.q, param_fit.pa, param_fit.x, param_fit.y,
            limit=-1)

    rad = pData['rad'] / param_fit.Re
    sbData = -2.5*log10(pData['sb']/0.396**2) + 22.5
    sbModel = -2.5*log10(pModel['sb']/0.396**2) + 22.5
    sbSigma = np.abs(2.5*sqrt(pVar['totalflux']) / pData['totalflux'] / log(10.) * sqrt(pVar['area']*0.396**2))

    ax4.plot(rad, sbData, 'o-', color='0.5')
    # ax4.errorbar(rad, sbData, yerr=sbSigma, fmt='ko-')
    good = sbSigma < 2  # to manage dynamical range
    ax4.fill_between(rad, sbData-sbSigma, sbData+sbSigma, where=good,
                     edgecolor='None', facecolor='0.85')
    ax4.plot(rad, sbModel, 'r-', lw=2)
    # ax4.axvline(rad[good][-1], ls='--', c='k')  # radius at which sigma(mag/arsec^2) reaches 2
    # ax_res.errorbar(rad, sbData-sbModel, yerr=sbSigma, fmt='bo-')
    ax_res.errorbar(rad, (sbData-sbModel)/sbSigma, fmt='bo-')
    ax_res.axhline(0, ls=':', c='k')
    # ax_res.plot(rad, pRes['mnflux'], 'ko-')
    # ax_res.plot(rad, pData['mnflux']-pModel['mnflux'], 'r-')
    mpltools.hide_tick_labels(ax4, 'x')
    ax4.invert_yaxis()
    ax4.minorticks_on()
    ax4.set_ylabel(r"$\mu$ (mag/arcsec$^2$)")
    ax_res.minorticks_on()
    ax_res.set_xlabel("$R/R_e$")
    ax_res.set_ylim(-5, 5)

    ax_chi = twinx(ax4)  #.twinx()
    radchi = pRedChi['rad']/param_fit.Re
    ax_chi.plot(radchi, pRedChi['mnflux'], 's-', c='b', mec='b', mfc='None')
    # cumulative reduced chi-squre
    cumChi = cumsum(pRedChi['totalflux'])
    cumN = cumsum(pRedChi['area'])
    cumRedChi = cumChi / cumN
    ax_chi.plot(radchi, cumRedChi, '.-', c='b')
    ax_chi.set_yscale('log')
    ax_chi.set_xlim(0, 8)
    ax_chi.axhline(1., ls='--', c='b')

    
    subplots_adjust(left=0.03, right=0.95, top=0.95, bottom=0.05, wspace=0.01)
    mpltools.hide_tick_labels(ax2, 'y')
    mpltools.hide_tick_labels(ax3, 'y')

    if 0:
        for ax in [ax1, ax2, ax3]:
            ax.set_xlim(param_fit.xc - param_fit.Re*5, param_fit.xc + param_fit.Re*5)
            ax.set_ylim(param_fit.yc - param_fit.Re*5, param_fit.yc + param_fit.Re*5)

    if savename:
        fig.savefig(savename)
        close()

    return vars()


def test_plot_residual():
    if 0:
        result = Table.read('sdss_psf/ser/deblended/RAWFIT00000.00048.fits')
        datadir = 'sdss_psf/ser/data/'
        modeldir = 'sdss_psf/ser/deblended/models'
        gal = result[9]
        v = plot_residual(gal['NAME'], Sersic(gal['FIT_SER']), datadir, modeldir)
    if 0:
        result = Table.read('sdss_psf/ser/deblended/RAWFIT00000.00048.fits')
        datadir = 'sdss_psf/ser/data/'
        modeldir = 'sdss_psf/ser/deblended/models'
        for i, gal in enumerate(result):
            print i
            v = plot_residual(gal['NAME'], Sersic(gal['FIT_SER']), datadir, modeldir,
                            savename='sdss_psf/ser/deblended/plotProfile/%s.png' % gal['NAME'])
    if 1:
        result = Table.read('sdss_psf/dvc/deblended/RAWFIT00000.00048.fits')
        datadir = 'sdss_psf/dvc/data/'
        modeldir = 'sdss_psf/dvc/deblended/models'
        for i, gal in enumerate(result):
            print i
            v = plot_residual(gal['NAME'], Sersic(gal['FIT_DVC']), datadir, modeldir,
                            savename='sdss_psf/dvc/deblended/plotProfile/%s.png' % gal['NAME'])
    return v


def main():
    """
    Save image/model/residual plots
    """
    import argparse
    parser = argparse.ArgumentParser(
                description="""Save image/model/residual plots""")
    parser.add_argument("input", type=str, help="part of NSA catalog")
    parser.add_argument("result", type=str, help="fit result table RAWXXX")
    parser.add_argument("profile", type=str, help="name of profile")
    parser.add_argument("datadir", type=str, help="input image directory")
    parser.add_argument("modeldir", type=str, help="model image directory")
    parser.add_argument("outdir", type=str, help="output directory")
    parser.add_argument("-r", "--range", nargs=2, type=int,
                        help="zero-based range of row index, e.g., -r 0 5")
    args = parser.parse_args()

    mpltools.mkdir_p(args.outdir)

    inputt = Table.read(args.input)
    result = Table.read(args.result)
    datadir, modeldir, profile, outdir = \
         args.datadir, args.modeldir, args.profile, args.outdir
    row_start = args.range[0] if args.range else 0
    row_end = args.range[1] + 1 if args.range else len(result)

    for index in range(row_start, row_end):
        gal = result[index]
        param_fit = Sersic(gal['FIT_%s' % (profile)])
        plot_residual(gal['NAME'], param_fit, datadir, modeldir,
            savename=outdir+'/%s.png'%(gal['NAME']))


if __name__ == '__main__':
    main()