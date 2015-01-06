""" Plotting fitting results """

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.io import fits
from asinh_norm import AsinhNorm
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from sersic import Sersic, NSersic, sersic2d, nanomaggie2mag
from .tools import radialprofile
from .fp import cosmo
import mpltools
from itertools import cycle

plt.rc('lines', lw=2)

__all__ = ["add_ellipse", "showim", "Model", "Sample"]


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
    e = Ellipse(xy=(x, y), width=2*rx, height=2*ry, angle=np.rad2deg(pa),
                **kwargs)
    ax.add_artist(e)


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
    ax_cb = divider.append_axes("bottom", size="5%", pad=0.2)
    cf = ax.get_figure()
    cf.add_axes(ax_cb)

    if not isinstance(norm, plt.Normalize):
        if norm is 'asinh':
            norm = AsinhNorm(
                # vmin=image[~isnan(image)].min(),
                vmin=np.percentile(image[~np.isnan(image)].flatten(), 5),
                vmax=np.percentile(image[~np.isnan(image)].flatten(), 90)*20)
        if norm is 'log':
            norm = LogNorm()
        if norm is None:
            pass
    im = ax.imshow(image, 
                   norm=norm,
                   interpolation='nearest',
                   aspect='auto', origin='lower', **kwargs)
    # ax.xaxis.tick_top()
    try:
        cb = plt.colorbar(im, cax=ax_cb, orientation='horizontal')
        # ax_cb.xaxis.tick_top()
        return ax, cb
    except:
        pass


class Model(object):
    """ Model class """
    def __init__(self, output_table_name, modeldir):
        self.output_table_name = output_table_name
        self.output_table = Table.read(output_table_name)
        self.modeldir = modeldir

        self.sky = self.output_table['SKY_%s' % (self.profile)]
        self.redchi = self.output_table['CHISQ_%s' % (self.profile)].data
        self.dof = self.output_table['DOF_%s' % (self.profile)].data
        self.chisq = self.redchi * self.dof
        self.p = NSersic(self.output_table['FIT_%s' % (self.profile)].data.T)

    @property
    def profile(self):
        self.output_table.colnames
        """ Get profile name from table column name """
        fitcol = [k for i, k in enumerate(self.output_table.colnames) if 'FIT_' in k]
        assert len(fitcol) == 1, "More than one profile found!"
        return fitcol[0].split('FIT_')[1]

    def get_image(self, ind):
        if 'OUTNAME' in self.output_table.colnames:
            imgname = self.modeldir+'/%s' % (self.output_table['OUTNAME'][ind])
        else:
            imgname = self.modeldir+'/M%s.fits' % (self.output_table['NAME'][ind])
        return fits.getdata(imgname)

    # @property
    # def ncomps(self):
    #     self.output_table[]


class Sample(object):
    """ sample class """
    def __init__(self, table_name, datadir, imgdir='images'):
        self.table_name = table_name
        self.table = Table.read(table_name)
        self.datadir = datadir
        self.models = {}
        self.imgdir = imgdir

    def add_model(self, name, table_name, modeldir):
        self.models[name] = Model(table_name, modeldir)

    def del_model(self, name):
        if name not in self.models.keys():
            print "cannot find %s" % (name)
        else:
            print "model %s removed" % (name)
            del self.models[name]

    def __str__(self):
        str = [
            'table : %s' % (self.table_name),
            'datadir : %s' % (self.datadir),
            'imgdir : %s' % (self.imgdir)
            ] + [self.models.__str__()]
        return '\n'.join(str)


    @property
    def Nmodels(self):
        return len(self.models.keys())

    def get_image(self, iden):
        return np.ma.masked_equal(fits.getdata(self.datadir + '/%s/%s.fits' % (self.imgdir, self.table['IAUNAME'][iden])), 0)

    def get_ivar(self, iden):
        return fits.getdata(self.datadir + '/ivar/%s.fits' % (self.table['IAUNAME'][iden]))

    def show_image(self, iden, unit='raw'):
        # calculate figure size
        image_size = 3.5  # inch
        margin = 0.2
        n_images = len(self.models.keys()) + 1  # models + data
        fig_x = n_images * (margin + image_size) + margin
        fig_y = image_size + 2*margin
        fig = plt.figure(figsize=(fig_x, fig_y))

        grid = plt.GridSpec(1, self.Nmodels+1)
        grid.update(left=0.05, right=0.97, top=0.9, bottom=0.1, wspace=0.02)
        # data image
        imgData = self.get_image(iden)
        imgIvar = self.get_ivar(iden)
        if unit is 'raw':
            ax, cax = showim(plt.subplot(grid[0]), imgData, norm='log', cmap='cubehelix')
        elif unit is 'mag':
            ax, cax = showim(plt.subplot(grid[0]), nanomaggie2mag(imgData/.396**2), cmap='cubehelix')
        # residual images
        for i, key in enumerate(self.models.keys(), start=1):
            imgModel = self.models[key].get_image(iden) + self.models[key].sky[iden]
            if unit is 'raw':
                ax, cax = showim(plt.subplot(grid[i]), (imgData - imgModel)*np.sqrt(imgIvar), cmap='RdBu', vmin=-5, vmax=5)
            elif unit is 'mag':
                resmag = nanomaggie2mag(imgData/.396**2) - nanomaggie2mag(imgModel/.396**2)
                ax, cax = showim(plt.subplot(grid[i]), resmag, cmap='cubehelix')
            ax.set_title(key)
            # highlight effective radius of each component
            p = self.models[key].p
            colors = cycle(['k', '#3CF5C1', '#F53C70'])
            for i in range(p.nprofiles):
                add_ellipse(ax, p[i].x[iden], p[i].y[iden], p[i].Re[iden],
                            p[i].Re[iden]*p[i].q[iden], p[i].pa[iden],
                            ec=colors.next())

    def get_profile_comps(self, iden, model):
        model = self.models[model]
        imgData = self.get_image(iden)
        center = self.table['CENTER'][iden]  # fixed galaxy center
        sersic = NSersic(model.output_table['FIT_%s'%(model.profile)][iden])
        x, y = np.meshgrid(np.arange(imgData.shape[1]),
                           np.arange(imgData.shape[0]))
        profs = []
        for i in range(sersic.nprofiles):
            r, p = radialprofile(sersic2d(x, y, sersic[i].p), center=center)
            profs.append((r, p))
        return profs

    def show_profile(self, iden, log=True, unit='raw', xunit='pixel', residual='sigma'):
        """
        show profile plot of a galaxy

        log : bool
            plot radius in log scale
        unit : str
            'raw' or 'mag'
        residual : str
            'sigma' or 'ratio'
        """
        assert residual in ['sigma', 'ratio'], "residual not recognized"
        assert unit in ['raw', 'mag'], "unit not recognized"
        assert xunit in ['pixel', 'kpc', 'arcsec'], "xunit not recognized"
        # setup figure and axes
        fig = plt.figure(figsize=(6, 6))
        fig.suptitle('%i %s' % (iden, self.table['IAUNAME'][iden][:6]), x=0.55, fontsize=15)
        fig.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.94, hspace=0.08)
        grid = plt.GridSpec(nrows=2, ncols=1)
        ax1 = plt.subplot(grid[0])
        ax2 = plt.subplot(grid[1], sharex=ax1)
        # color cycle
        colors = cycle(['#348ABD', '#A60628', '#7A68A6', '#467821',
                  '#D55E00', '#CC79A7', '#56B4E9', '#009E73', '#F0E442', '#0072B2'])

        # get data
        center = self.table['CENTER'][iden]  # fixed galaxy center
        imgData = self.get_image(iden)
        rData, pData = radialprofile(imgData, center=center, binsize=2, mask=imgData.mask)
        rData, sigData = radialprofile(imgData, statistic=np.std, center=center, binsize=2, mask=imgData.mask)
        rData, countData = radialprofile(imgData, statistic='count', center=center, binsize=2, mask=imgData.mask)
        
        # radius conversion factor
        if xunit is 'pixel':
            xconv = 1.
            ax2.set_xlabel("radius (pixels)")
        elif xunit is 'arcsec':
            xconv = 0.396
            ax2.set_xlabel("radius (arcsec)")
        elif xunit is 'kpc':
            xconv = 0.396 * cosmo.angular_diameter_distance(self.table['Z_1'][iden]).value * 1000./206265.
            ax2.set_xlabel("radius (kpc)")
        # start plotting
        if unit is 'raw':
            ax1.plot(rData*xconv, pData, 'k.', zorder=5)
            ax1.set_ylabel("nanomaggies / pixel")
        elif unit is 'mag':
            ax1.plot(rData*xconv, nanomaggie2mag(pData/.396**2), 'k.', zorder=5)
            ax1.set_ylabel("mag / arcsec$^2$")
        # plot each model
        jmodel = 1
        for key, model in self.models.iteritems():
            c = colors.next()
            rModel, pModel = radialprofile(model.get_image(iden) + model.sky[iden], center=center, binsize=2)
            if unit is 'raw':
                ax1.plot(rModel*xconv, pModel, label=key, c=c)
            elif unit is 'mag':
                ax1.plot(rModel*xconv, nanomaggie2mag(pModel/.396**2), label=key, c=c)
            if residual is 'sigma':
                ax2.set_ylabel('residual / sigma')
                ax2.plot(rModel*xconv, (pData-pModel)/(sigData/np.sqrt(countData)), c=c)
                ax2.axhline(0, c='k', lw=1)
            elif residual is 'ratio':
                ax2.plot(rModel*xconv, pData/pModel, c=c)
            # indicate sub-components
            sersic = NSersic(model.output_table['FIT_%s'%(model.profile)][iden])
            if sersic.nprofiles > 1:
                x, y = np.meshgrid(np.arange(imgData.shape[1]),
                                   np.arange(imgData.shape[0]))
                for i in range(sersic.nprofiles):
                    r, p = radialprofile(sersic2d(x, y, sersic[i].p), center=center)
                    if unit is 'raw':
                        ax1.plot(r*xconv, p, ls='--', lw=1, c=c)
                    elif unit is 'mag':
                        ax1.plot(r*xconv, nanomaggie2mag(p/.396**2), ls='--', lw=1, c=c)
            # indicate effective radius
            for i in range(sersic.nprofiles):
                ax1.axvline(sersic[i].Re*xconv, color=c, lw=1, ymin=0.05, ymax=0.05+0.05*jmodel)
            jmodel += 1

        ax1.set_xlim([.95*rData.min()*xconv, 1.05*rData.max()*xconv])
        plt.setp(ax1.get_xticklabels(), visible=False)
        if unit is 'raw':
            ax1.set_yscale('log')
            ax1.set_ylim(pData[pData>0].min(), pData.max())
        elif unit is 'mag':
            ax1.set_ylim(nanomaggie2mag(pData[pData>0].min()*0.9/.396**2), nanomaggie2mag(pData.max()*1.1/.396**2))
        ax1.legend(loc='upper right')
        # ax2.set_xlim([0, xlim/2.])

        if log is True:
            ax1.set_xscale('log')
            ax2.set_xscale('log')



