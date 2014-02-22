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

class Parameters(object):
    """ single Sersic profile parameters"""
    def __init__(self, param):
       self.Ie = param[0]
       self.Re = param[1]
       self.N = param[2]
       self.ab = param[3]
       self.shape = param[4]
       self.xc = param[5]
       self.yc = param[6]
       self.posang = param[7]

outdir = './plots'
datadir = './data'
modeldir = './out_SER/models'
inputt = Table.read('input_SER.fits')
result = Table.read('./out_SER/RAWFIT00000.00050.fits')

# for index in arange(20, 30):
for index in [20]:
    profile = 'SER'

    # load images
    imagefits = datadir+'/images/%s-%i-atlas-%i.fits.gz' % (
            result[index]['NAME'], result[index]['PARENT_ID'], result[index]['ATLAS_ID'])
    modelfits = modeldir+'/M%s.fits' % (result[index]['NAME'])
    model = fits.getdata(modelfits, ext=0)
    data  = fits.getdata(imagefits, ext=2)  # r-band
    residual = data - model

    if 0:
        fig = figure(figsize=(10, 4))
        fig.suptitle('%3i %s' % (index, result[index]['NAME']))
        ax1 = fig.add_subplot(1,3,1)
        ax1.set_title('data')
        ax1.imshow(data)
        ax2 = fig.add_subplot(1,3,2)
        ax2.set_title('model')
        ax2.imshow(model)
        fig.add_subplot(1,3,3)
        title('residual')
        imshow(residual)
        # colorbar()
        tight_layout()

    iauname = result[index]['NAME']
    param_init = Parameters(result[index]['%s_VAL' % profile])
    param_fit  = Parameters(result[index]['FIT_%s' % profile])
    xw0, yw0

    # use aplpy
    fig = figure(figsize=(12, 5))
    fig.suptitle('%3i %s' % (index, iauname))

    vmin = min([min(r) for r in data])
    vmax = max([max(r) for r in data])
    vmid = vmax/2.

    # data image
    f1 = aplpy.FITSFigure(imagefits, hdu=2, figure=fig, subplot=[0.05,0.2,0.3,0.7])
    f1.show_grayscale(stretch='arcsinh')
    f1.axis_labels.hide()
    f1.tick_labels.hide()
    # add initial and best-fit center
    f1.show_markers(*f1.pixel2world(param_fit.xc, param_fit.yc),
                    edgecolor='r', facecolor='r', lw=1)
    f1.show_markers(*f1.pixel2world(param_init.xc, param_init.yc),
                    edgecolor='c', facecolor='c', lw=1)
    # add initial effective radius ellipse
    f1.show_ellipses(
            param_init.xc, param_init.yc, param_init.Re*2, param_init.Re*2./param_init.ab,
            angle=rad2deg(param_init.posang))

    # best-fit model
    f2 = aplpy.FITSFigure(modelfits, hdu=0, figure=fig, subplot=[0.36,0.2,0.3,0.7])
    f2.show_grayscale(stretch='arcsinh')
    f2.axis_labels.hide()
    f2.tick_labels.hide()
    # add initial and best-fit center
    f2.show_markers(*f2.pixel2world(param_fit.xc, param_fit.yc),
                    edgecolor='r', facecolor='r', lw=1)
    f2.show_markers(*f2.pixel2world(param_init.xc, param_init.yc),
                    edgecolor='c', facecolor='c', lw=1)
    # add model effective radius ellipse
    f2.show_ellipses(
            param_fit.xc, param_fit.yc, param_fit.Re*2, param_fit.Re*2*param_fit.ab,
            angle=rad2deg(param_fit.posang))

    # residual image
    f3 = aplpy.FITSFigure(residual, figure=fig, subplot=[0.67,0.2,0.3,0.7])
    f3.show_grayscale()
    f3.axis_labels.hide()
    f3.tick_labels.hide()

    # info = []
    # for i, p in enumerate(result[index]['FIT_%s' % profile]):
    #     fixed = '*' if result[index]['%s_FIX' % profile][i] else ''
    #     info.append(fixed + '%s = %f' % (pname[i], p))
    # info = []
    # info_str = ' '.join(info[:4]) + '\n' + ' '.join(info[4:])
    # figtext(0.5, 0.08, info_str, va='center', ha='center')

    show()
    raw_input()
    

