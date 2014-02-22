""" Script to plot fitting results """

from pylab import *
from astropy.table import Table
from astropy.io import fits
import aplpy

def sersic_param(param):
    if len(param) % 8 != 0:
        raise 'param array needs to be a multiple of 8'
    param = atleast_1d(param)
    Nprofiles = param.size / 8
    dic = {}
    for i, p in enumerate(param.reshape([Nprofiles, 8])):
        print '[%i]' % i, p

# parameter dicts
pind = {
    'I0': 0,
    'Reff': 1,
    'index': 2,
    'ratio': 3,
    'shape': 4,
    'center_x': 5,
    'center_y': 6,
    'posang': 7,
}
pname = {v:k for k, v in pind.items()}


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


class Parameters(AttrDict):
    """ single Sersic profile parameters"""
    def __init__(self, param):
        d = {}
        d['Ie'] = param[0]
        d['Re'] = param[1]
        d['N'] = param[2]
        d['ba'] = param[3]
        d['shape'] = param[4]
        d['xc'] = param[5]
        d['yc'] = param[6]
        d['posang'] = param[7]
        AttrDict.__init__(self, d)

    def __str__(self):
        info = []
        info.append('%12s = %10.6f' % ('Ie', self.Ie))
        info.append('%12s = %10.6f' % ('Re', self.Re))
        info.append('%12s = %10.6f' % ('N', self.N))
        info.append('%12s = %10.6f' % ('ba', self.ba))
        info.append('%12s = %10.6f' % ('xc', self.xc))
        info.append('%12s = %10.6f' % ('yc', self.yc))
        info.append('%12s = %10.6f' % ('posang', self.posang))
        return '\n'.join(info)

    def __repr__(self):
        return self.__str__()


def showim(ax, image):
    from asinh_norm import AsinhNorm

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    ax_cb = divider.append_axes("bottom", size="5%", pad=0.05)
    fig1 = ax.get_figure()
    fig1.add_axes(ax_cb)

    vmin = min([min(r) for r in image])
    im = ax.imshow(image, cmap=cm.gray_r, 
            norm=matplotlib.colors.LogNorm(),
            #interpolation='nearest',
            aspect='auto')
    ax.xaxis.tick_top()
    try:
        cb = colorbar(im, cax=ax_cb, orientation='horizontal')
        ax_cb.xaxis.tick_bottom()
        return ax, cb
    except:
        pass

    


outdir = './plots'
datadir = './data'
modeldir = './out_SER2/models'
inputt = Table.read('SampleZMprobaEllSub_visual.fits')
result = Table.read('./out_SER2/RAWFIT00000.00004.fits')

# for index in arange(20, 30):
for index in arange(0, 5):
    profile = 'SER'

    # load images
    imagefits = datadir+'/images/%s-%i-atlas-%i.fits.gz' % (
            result[index]['NAME'], result[index]['PARENT_ID'], result[index]['ATLAS_ID'])
    modelfits = modeldir+'/M%s.fits' % (result[index]['NAME'])
    model = fits.getdata(modelfits, ext=0)
    data  = fits.getdata(imagefits, ext=2)  # r-band
    residual = data - model

    iauname = result[index]['NAME']
    z = inputt[index]['Z_1']
    param_init = Parameters(result[index]['%s_VAL' % profile])
    param_fit  = Parameters(result[index]['FIT_%s' % profile])

    # crop images
    # halfsize = param_init.Re * 3
    # data = data[param_init.yc-halfsize:param_init.yc+halfsize,
    #         param_init.xc-halfsize:param_init.xc+halfsize]
    # model = model[param_init.yc-halfsize:param_init.yc+halfsize,
    #                 param_init.xc-halfsize:param_init.xc+halfsize]
    # residual = residual[param_init.yc-halfsize:param_init.yc+halfsize,
    #                 param_init.xc-halfsize:param_init.xc+halfsize]

    # use aplpy
    fig = figure(figsize=(14, 5))
    fig.suptitle('%3i %s %6.4f' % (index, iauname, z))

    transData2Axes = lambda ax, p: ax.transAxes.inverted().transform(ax.transData.transform(p))

    vmin = min([min(r) for r in data])
    vmax = max([max(r) for r in data])
    vmid = vmax/2.

    # data image
    ax1 = fig.add_subplot(1,3,1)
    ax1.text(0.5, 1.1, 'data', ha='center', va='bottom', transform=ax1.transAxes)
    showim(ax1, data)
    # add initial and best-fit center
    ax1.axvline(param_init.xc, c='c')
    ax1.axhline(param_init.yc, c='c')
    ax1.axvline(param_fit.xc, c='m')
    ax1.axhline(param_fit.yc, c='m')
    # add initial effective radius ellipse
    from matplotlib.patches import Ellipse
    e = Ellipse(xy=(param_init.xc, param_init.yc), width=2*param_init.Re,
                    height=2*param_init.Re/param_init.ba, angle=rad2deg(param_init.posang)+90.,
                    fc='None', ec='c')
    e2 = Ellipse(xy=(param_fit.xc, param_fit.yc), width=2*param_fit.Re,
                    height=2*param_fit.Re/param_fit.ba, angle=rad2deg(param_fit.posang)+90.,
                    fc='None', ec='m')
    ax1.add_artist(e)
    ax1.add_artist(e2)

    # best-fit model
    ax2 = fig.add_subplot(1,3,2, sharex=ax1, sharey=ax1)
    ax2.text(0.5, 1.1, 'model', ha='center', va='bottom', transform=ax2.transAxes)
    showim(ax2, model)
    # add initial and best-fit center
    ax2.axvline(param_fit.xc, c='m')
    ax2.axhline(param_fit.yc, c='m')

    # residual image
    ax3 = fig.add_subplot(1,3,3, sharex=ax1, sharey=ax1)
    ax3.text(0.5, 1.1, 'residual', ha='center', va='bottom', transform=ax3.transAxes)
    showim(ax3, residual)

    
    figtext(0.95, 0.81, '%r' % param_fit, va='top', ha='right',
            family='monospace')

    tight_layout(rect=(0.01,0.1,0.814,0.9))
    # show()
    fig.savefig('out_SER2/%s.png' % iauname)
    

