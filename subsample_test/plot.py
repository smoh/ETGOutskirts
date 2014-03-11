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
    parser.add_argument("datadir", type=str, help="data directory")
    parser.add_argument("modeldir", type=str, help="model image directory")
    parser.add_argument("outdir", type=str, help="output directory")
    parser.add_argument("-r", "--range", nargs='+', type=int,
                        help="zero-based range of row index, e.g., -r 0 5")
    args = parser.parse_args()

    inputt = Table.read(args.input)
    result = Table.read(args.result)
    datadir, modeldir, profile, outdir = args.datadir, args.modeldir, \
            args.profile, args.outdir
    row_start = args.range[0] if args.range else 0
    row_end = args.range[1] if args.range else len(inputt) + 1

    for index in range(row_start, row_end):

        # load images
        imagefits = datadir+'/images/%s-%i-atlas-%i.fits.gz' % (
                result[index]['NAME'], result[index]['PARENT_ID'], result[index]['ATLAS_ID'])
        modelfits = modeldir+'/M%s.fits' % (result[index]['NAME'])
        model = fits.getdata(modelfits, ext=0)
        data  = fits.getdata(imagefits, ext=2)  # r-band
        residual = data - model

        # galaxy info
        nsaid = inputt[index]['NSAID']
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

        fig = figure(figsize=(14, 5))
        fig.suptitle('%3i NSAID %6i %s %6.4f' % (index, nsaid, iauname, z))

        transData2Axes = lambda ax, p: ax.transAxes.inverted().transform(ax.transData.transform(p))

        vmin = min([min(r) for r in data])
        vmax = max([max(r) for r in data])
        vmid = vmax/2.

        def add_effective_ellipse(ax):
            from matplotlib.patches import Ellipse
            e = Ellipse(xy=(param_init.xc, param_init.yc), width=2*param_init.Re,
                            height=2*param_init.Re/param_init.ba, angle=rad2deg(param_init.posang)+90.,
                            fc='None', ec='c')
            e2 = Ellipse(xy=(param_fit.xc, param_fit.yc), width=2*param_fit.Re,
                            height=2*param_fit.Re/param_fit.ba, angle=rad2deg(param_fit.posang)+90.,
                            fc='None', ec='c')
            # ax.add_artist(e)
            ax.add_artist(e2)

        # data image
        ax1 = fig.add_subplot(1,3,1)
        ax1.text(0.5, 1.1, 'data', ha='center', va='bottom', transform=ax1.transAxes)
        showim(ax1, data)
        if 0:
            # add initial and best-fit center crosshair
            ax1.axvline(param_init.xc, c='c')
            ax1.axhline(param_init.yc, c='c')
            ax1.axvline(param_fit.xc, c='m')
            ax1.axhline(param_fit.yc, c='m')
        add_effective_ellipse(ax1)

        # best-fit model
        ax2 = fig.add_subplot(1,3,2, sharex=ax1, sharey=ax1)
        ax2 = fig.add_subplot(1,3,2, sharex=ax1, sharey=ax1)
        ax2.text(0.5, 1.1, 'model', ha='center', va='bottom', transform=ax2.transAxes)
        showim(ax2, model)
        add_effective_ellipse(ax2)

        # residual image
        ax3 = fig.add_subplot(1,3,3, sharex=ax1, sharey=ax1)
        ax3.text(0.5, 1.1, 'residual', ha='center', va='bottom', transform=ax3.transAxes)
        showim(ax3, residual)
        add_effective_ellipse(ax3)

        
        # figtext(0.95, 0.81, '%r' % param_fit, va='top', ha='right',
        #         family='monospace')

        subplots_adjust(left=0.05, right=0.95, top=0.80, bottom=0.1)
        # tight_layout(rect=(0.01,0.1,0.814,0.9))
        fig.savefig(outdir + '/%s.png' % (iauname))
        close()


if __name__ == '__main__':
    main()
