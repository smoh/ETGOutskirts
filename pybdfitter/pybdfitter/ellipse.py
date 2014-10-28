"""
collection of tools to interact with iraf ellipse
"""

from pyraf import iraf
import os, shutil, uuid, re

__all__ = ['ellipse', 'read_ellipse_output']

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
    # load packages
    iraf.stsdas(_doprint=0)
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
        filename, format='ascii.commented_header', header_start=1, fill_values=['INDEF', '-99'])
