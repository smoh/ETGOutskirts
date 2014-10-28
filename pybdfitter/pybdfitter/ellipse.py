"""
collection of tools to interact with iraf ellipse
"""

import os, shutil, uuid, re
from pyraf import iraf
from astropy.table import Table
from matplotlib import pyplot as plt

__all__ = ['ellipse', 'read_ellipse_output', 'EllipsePlot']

def shorten_iauname(name):
    p = re.compile('(J[\d.]+)[+-][\d.]+')
    matches = p.findall(name)
    return matches[0][:6]

def ellipse(image, outname, verbose=False, **kwargs):
    """run ellipse task

    image : input image name
    outname : output filename
    
    Keywords:
    *geompar*
        x0, y0 -- required
        ellip0 : initial ellipticity (default=0.2)
        pa0 : initial PA in degrees ccw from +y (default=20)
        sma0 : initial semi-major axis (default=10)
        minsma : minimum sma (default=1)
        maxsma : maximum sma (default='INDEF')
        step : sma step size (default=0.1)
        linear : linear sampling of sma? (default='no')
    *controlpar*
        conver : 
    """
    # verify input
    if 'x0' not in kwargs.keys() or 'y0' not in kwargs.keys():
        raise KeyError('x0 or y0 is not found')
    # load packages
    iraf.stsdas(_doprint=0, motd=False)
    iraf.analysis(_doprint=0)
    iraf.isophote(_doprint=0)

    # reset all parameters
    iraf.unlearn('ellipse', 'geompar', 'controlpar', 'samplepar')
    iraf.ellipse.interactive = False

    # geompar
    iraf.ellipse.x0 = kwargs.pop('x0')
    iraf.ellipse.y0 = kwargs.pop('y0')
    iraf.ellipse.ellip0 = kwargs.pop('ellip0', 0.2)
    iraf.ellipse.pa0 = kwargs.pop('pa0', 20.0)
    iraf.ellipse.sma0 = kwargs.pop('sma0', 10.0)
    iraf.ellipse.minsma = kwargs.pop('minsma', 1.)
    iraf.ellipse.maxsma = kwargs.pop('maxsma', 'INDEF')
    iraf.ellipse.step = kwargs.pop('step', 0.1)
    iraf.ellipse.linear = kwargs.pop('linear', 0)
    iraf.ellipse.verbose = verbose

    # controlpar
    iraf.ellipse.conver = kwargs.pop('conver', 0.05)
    iraf.ellipse.minit = kwargs.pop('minit', 10)
    iraf.ellipse.maxit = kwargs.pop('maxit', 50)
    iraf.ellipse.hcenter = kwargs.pop('hcenter', 'no')

    # samplepar
    iraf.ellipse.tsample = kwargs.pop('tsample', 'none')
    iraf.ellipse.harmonics = kwargs.pop('harmonics', 'none')

    if verbose:
        iraf.lparam('ellipse')
        iraf.lparam('geompar')
        iraf.lparam('controlpar')
        iraf.lparam('samplepar')

    # run ellipse
    tempname = shorten_iauname(image) + str(uuid.uuid4()) + '.tab'
    iraf.ellipse(image, tempname)
    # print ascii table
    iraf.tprint(tempname, pwidth='INDEF', Stdout=outname)
    os.remove(tempname)

def read_ellipse_output(filename):
    return Table.read(
        filename, format='ascii.commented_header', header_start=1, data_start=0,
        fill_values=['INDEF', '-99'])


class EllipsePlot(object):
    def __init__(self, out):
        self.fn = out
        self.t = read_ellipse_output(out)
        self.badexit = self.t['STOP'] != 0
        self.good = self.t['RMS'] < self.t['INTENS']

    def setup_figure(self):
        fig = plt.figure(figsize=(6, 8))
        fig.subplots_adjust(left=0.18, right=0.95, wspace=.4, top=.95)
        axIntens = plt.subplot2grid((4,1), (0,0), rowspan=2)
        axPA     = plt.subplot2grid((4,1), (2,0), sharex=axIntens)
        axEllip  = plt.subplot2grid((4,1), (3,0), sharex=axIntens)

        self.fig = fig
        self.plot_intens(axIntens)
        self.plot_pa(axPA)
        self.plot_ellip(axEllip)

    def plot_intens(self, ax=None):
        t = self.t
        ax.errorbar(t['SMA'], t['INTENS'], yerr=t['RMS'], fmt='o-')
        ax.plot(t['SMA'][self.badexit], t['INTENS'][self.badexit], 'ks', mfc='None', ms=12)
        ax.set_yscale('log')

    def plot_pa(self, ax=None):
        t = self.t
        pa = t['PA']
        # pa[pa<0] = pa[pa<0] +180.
        ax.errorbar(t['SMA'], pa, yerr=t['PA_ERR'], fmt='o-')
        ax.plot(t['SMA'][self.badexit], pa[self.badexit], 'ro', mfc='None', ms=12)
        ax.set_ylabel('PA')

    def plot_ellip(self, ax=None):
        t = self.t
        ax.errorbar(t['SMA'], t['ELLIP'], yerr=t['ELLIP_ERR'], fmt='o-')
        ax.plot(t['SMA'][self.badexit], t['ELLIP'][self.badexit], 'ro', mfc='None', ms=12)
        ax.plot(t['SMA'], t['A_BIG']/t['RMS'], 'ko:')
        ax.set_ylabel('Ellipticity')

    def plot_3(self, ax=None):
        t = self.t
        ax.errorbar(t['SMA'], t['A3'], yerr=t['A3_ERR'], fmt='o-')
        ax.errorbar(t['SMA'], t['B3'], yerr=t['B3_ERR'], fmt='o-')
        ax.set_ylabel('A3, B3')

    def plot_4(self, ax=None):
        t = self.t
        ax.errorbar(t['SMA'], t['A4'], yerr=t['A4_ERR'], fmt='o-')
        ax.errorbar(t['SMA'], t['B4'], yerr=t['B4_ERR'], fmt='o-')
        ax.set_ylabel('A4, B4')
