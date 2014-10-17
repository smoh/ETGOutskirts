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
    im = ax.imshow(image, cmap=plt.cm.jet,
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
        return fits.getdata(self.modeldir+'/M%s.fits' % (self.output_table['NAME'][ind]))


    # @property
    # def ncomps(self):
    #     self.output_table[]


class Sample(object):
    """ sample class """
    def __init__(self, table_name, datadir):
        self.table_name = table_name
        self.table = Table.read(table_name)
        self.datadir = datadir
        self.models = {}

    def add_model(self, name, table_name, modeldir):
        self.models[name] = Model(table_name, modeldir)

    def del_model(self, name):
        if name not in self.models.keys():
            print "cannot find %s" % (name)
        else:
            print "model %s removed" % (name)
            del self.models[name]

    @property
    def Nmodels(self):
        return len(self.models.keys())

    def get_image(self, iden):
        return np.ma.masked_equal(fits.getdata(self.datadir + '/images_admask/%s.fits' % (self.table['IAUNAME'][iden])), 0)

    def get_ivar(self, iden):
        return fits.getdata(self.datadir + '/ivar/%s.fits' % (self.table['IAUNAME'][iden]))

    def show_image(self, iden):
        fig = plt.figure(figsize=(15, 3.5))
        grid = plt.GridSpec(1, self.Nmodels+1)
        grid.update(left=0.05, right=0.97, top=0.9, bottom=0.1, wspace=0.02)
        # data image
        imgData = self.get_image(iden)
        ax, cax = showim(plt.subplot(grid[0]), imgData, norm='log')
        # residual images
        for i, key in enumerate(self.models.keys(), start=1):
            imgModel = self.models[key].get_image(iden) + self.models[key].sky[iden]
            ax, cax = showim(plt.subplot(grid[i]), imgData - imgModel, norm='log')
            ax.set_title(key)
            # highlight effective radius of each component
            p = self.models[key].p
            colors = cycle(['k', '#3CF5C1', '#F53C70'])
            for i in range(p.nprofiles):
                add_ellipse(ax, p[i].x[iden], p[i].y[iden], p[i].Re[iden],
                            p[i].Re[iden]*p[i].q[iden], p[i].pa[iden],
                            ec=colors.next())

    def show_profile(self, iden):
        # setup figure and axes
        fig = plt.figure(figsize=(6, 8))
        fig.suptitle('%i %s' % (iden, self.table['IAUNAME'][iden][:6]), x=0.55, fontsize=15)
        fig.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.94, hspace=0.08)
        grid = plt.GridSpec(nrows=2, ncols=1)
        ax1 = plt.subplot(grid[0])
        ax2 = plt.subplot(grid[1], sharex=ax1)
        center = self.table['CENTER'][iden]  # fixed galaxy center
        imgData = self.get_image(iden)
        rData, pData = radialprofile(imgData, center=center, binsize=2, mask=imgData.mask)
        rData, sigData = radialprofile(imgData, statistic=np.std, center=center, binsize=2, mask=imgData.mask)
        rData, countData = radialprofile(imgData, statistic='count', center=center, binsize=2, mask=imgData.mask)
        ax1.plot(rData, pData, 'ko')
        colors = cycle(['#348ABD', '#A60628', '#7A68A6', '#467821',
                  '#D55E00', '#CC79A7', '#56B4E9', '#009E73', '#F0E442', '#0072B2'])
        residual = 'sigma'  # which residual plot?
        # plot each model
        for key, model in self.models.iteritems():
            c = colors.next()
            rModel, pModel = radialprofile(model.get_image(iden) + model.sky[iden], center=center, binsize=2)
            ax1.plot(rModel, pModel, label=key, c=c)
            # ax2.plot(rModel, pModel/pData)
            if residual is 'sigma':
                ax2.plot(rModel, (pData-pModel)/(sigData/np.sqrt(countData)), c=c)
            elif residual is 'ratio':
                ax2.plot(rModel, pData/pModel, c=c)

            # indicate sub-components
            sersic = NSersic(model.output_table['FIT_%s'%(model.profile)][iden])
            if sersic.nprofiles > 1:
                x, y = np.meshgrid(np.arange(imgData.shape[1]),
                                   np.arange(imgData.shape[0]))
                for i in range(sersic.nprofiles):
                    r, p = radialprofile(sersic2d(x, y, sersic[i].p), center=center)
                    ax1.plot(r, p, ls='--', c=c)
                # clist = cycle(['g', 'm'])

        # ax1.set_xscale('log')
        xlim = ax1.get_xlim()[1]
        # ax1.set_xlim([0, xlim/2.])
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.set_yscale('log')
        ax1.set_ylim(pData[pData>0].min(), pData.max())
        ax1.set_ylabel("nanomaggies / pixel")
        ax1.legend(loc='upper right')

        # ax2.set_xlim([0, xlim/2.])
        ax2.set_xlabel("radius (pixels)")
        if residual is 'sigma':
            ax2.set_ylabel('residual / sigma')


        ax2.axhline(1, c='k')

