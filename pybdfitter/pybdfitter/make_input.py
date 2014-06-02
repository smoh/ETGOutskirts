#!/usr/bin/env python
"""
Script to make input FITS for bdfitter

-------------
Claire's note
-------------
[NAME OF PROFILE]_VAL
[NAME OF PROFILE]_FIX
each with 8 X (the number of profiles)

The name of the profile is whatever you want, but you need to give it
to fit_sample, so it can find the input parameters
XXX_FIX is a boolean array which =1 for fixed parameters
XXX_VAL is the input values, including fixed values

if you don't fix the initial position of flux, default values will 
be chosen by the code (based on image size)
the parameters (in order) are:

0 -- the surface brightness
1 -- the half light radius (refers to the semimajor axis)
2 -- Sersic index
3 -- axis ratio (minor/major)
4 -- shape of isophote, default is zero and fixed, 
        but if you specify XXX_FIX you need to fix the 4th parameter
5 -- x coordinate of center of profile, be default multiple profiles 
        share a center
6 -- y coordinate of center
7 -- position angle (radians counterclockwise from x-axis)

you need to pick starting values for the size, and the axis ratio
starting values for the surface brightness and the position will
be set in the code, don't set the flux to zero, as the code
just rescales values (so can't rescale 0)

for two component fits, always but the *bulge* (smaller +higher sersic)
profile first
you'll have to source the fixed values from somewhere, either
based on the radius, or based on previously fitting the profiles
                    
it might be good to have different input scripts, for example
you could make 1 input file to fit the sersic/dvc profile and after
that runs, make another file for the two component fits, which 
takes the olds fits as inputs for the VALs arrays
    
"""

from pylab import *
import sys, os, re
import argparse
import numpy as np
import pyfits

from astropy.table import Table, Column
from astropy.io import fits
import argparse

from pybdfitter.sersic import ip, nanomaggie2mag


class InputTable(Table):
    """Input table for bdfitter"""
    def __init__(self, source):
        super(InputTable, self).__init__()
        self.datadir = 'data'
        if type(source) == str:
            self.srctab = Table.read(source)
        else:
            self.srctab = source
        self._get_default_entries()

    def _get_default_entries(self):
        """
        IAUNAME : identifier
        RUN, CAMCOL, FIELD, XPOS, YPOX : for reading PSF
        SPA : for rotating PSF
        """
        if 'IAUNAME' in self.srctab.colnames:
            self.add_column(self.srctab['IAUNAME'])
            self['IAUNAME'].name = 'NAME'
        elif 'NAME' in self.srctab.colnames:
            self.add_column(self.srctab['NAME'])
        else:
            raise InputError, 'no IAUNAME or NAME column found'
        for k in ['RUN', 'CAMCOL', 'FIELD', 'XPOS', 'YPOS', 'VDISP', 'Z_1']:
            if not k in self.srctab.colnames:
                raise ValueError, "%s not found" % (k)
            self.add_column(self.srctab[k])

        # rotation angle for SDSS PSF
        angle = []
        if 'SPA' in self.srctab.colnames:
            self.add_column(self.srctab['SPA'])
        else:
            for row in self:
                name = self.datadir + '/sdss_field/frame-r-%06i-%i-%04i.fits' % (
                    row['RUN'], row['CAMCOL'], row['FIELD'])
                angle.append(fits.getheader(name)['SPA'])
            self.add_column(Column(name='SPA', data=array(angle).astype(float)))

    def add_profile(self, name, ncomp):
        """
        name : name of the profile
        ncomp : the number of components
        """
        Nrows = len(self.srctab)
        val = Column(
            name='%s_VAL' % (name),
            data=array([[1, 10, 4, 0.7, 0, 0, 0, 0]*ncomp]*Nrows).astype(float))
        fix = Column(
            name='%s_FIX' % (name),
            data=array([[0, 0, 0, 0, 1, 0, 0, 0]*ncomp]*Nrows))
        self.add_columns([val, fix])

        
    def set_values_from_NSA(self, name, icomp):
        """
        Set initial values of icomp-th component of profile to
        those given by NSA single Sersic profile

        name : name of the profile
        icomp : zero-based i-th component to set the values
        """
        # get initial values from NSA Sersic fitting
        val = '%s_VAL' % (name)
        self[val][:,icomp*8 + ip['pa']] = deg2rad((self.srctab['SERSIC_PHI'].data + 90.) % 360.)
        self[val][:,icomp*8 + ip['q']] = self.srctab['SERSIC_BA'].data
        self[val][:,icomp*8 + ip['Re']] = self.srctab['SERSIC_TH50'].data / 0.396
        self[val][:,icomp*8 + ip['x']] = self.srctab['XCEN'].data
        self[val][:,icomp*8 + ip['y']] = self.srctab['YCEN'].data     

    def write(self, *args, **kwargs):
        if not kwargs.has_key('overwrite'):
            kwargs['overwrite'] = True
        # ugly way to get registry working
        Table(self).write(*args, **kwargs)   



def single_SER(input, outname):
    """
    generate input table for single Sersic profile

    input : str
        input table name; should be part of NSA catalog in FITS format
    outname : str
        output table name
    """
    t = InputTable(input)
    t.add_profile('SER', 1)
    t.set_values_from_NSA('SER', 0)
    t.write(outname)


def single_DVC(input, outname):
    """
    generate input table for single de Vaucouleurs profile

    input : str
        input table name; should be part of NSA catalog in FITS format
    outname : str
        output table name
    """
    t = InputTable(input)
    t.add_profile('DVC', 1)
    t.set_values_from_NSA('DVC', 0)
    # fix sersic index
    t['DVC_VAL'][:,ip['n']] = 4.
    t['DVC_FIX'][:,ip['n']] = 1.
    t.write(outname, format='fits')


def dvc_exp(input, outname):
    """
    DVC + EXP profile
    """
    t = InputTable(input)
    t.add_profile('DVCEXP', 2)
    t.set_values_from_NSA('DVCEXP', 0)
    t.set_values_from_NSA('DVCEXP', 1)
    # fix Sersic indices
    t['DVCEXP_VAL'][:,ip['n']] = 4.
    t['DVCEXP_FIX'][:,ip['n']] = 1.
    t['DVCEXP_VAL'][:,8+ip['n']] = 1.
    t['DVCEXP_FIX'][:,8+ip['n']] = 1.
    t['DVCEXP_VAL'][:,8+ip['Re']] = 2*t['DVCEXP_VAL'][:,ip['Re']]
    t['DVCEXP_VAL'][:,8+ip['Ie']] = 0.1*t['DVCEXP_VAL'][:,ip['Ie']]
    t.write(outname, format='fits')




# TODO: functions below need to be modified
def single_DVC_fixRe(input, outname):
    """
    generate input table for single De Vauc profile, with
    circularized effective radius fixed

    inputtable : str
        input table name; should be part of NSA catalog in FITS format
    outname : str
        output table name
    """
    t = InputTable(input)
    t.add_profile('DVC', 1)
    t.set_values_from_NSA('DVC', 0)
    # fix sersic index
    t['DVC_VAL'][:,ip['n']] = 4.
    t['DVC_FIX'][:,ip['n']] = 1.

    # effective radius from FP
    from fp import Bernardi03_r
    sigma = t.srctab['VDISP']
    flux_mag = nanomaggie2mag(t.srctab['PETROFLUX'][:,4])
    z = t.srctab['Z_1']
    t['DVC_VAL'][:,ip['Re']] = Bernardi03_r.r0(sigma, flux_mag, z)/0.396 /t['DVC_VAL'][:,ip['q']]
    # t.write(outname, clobber=True)



def double_SER(outname):
    """
    generate input table for double Sersic profile

    inputtable : str
        input table name; should be part of NSA catalog in FITS format
    outname : str
        output table name
    """
    # this will be arguments at some point
    input = 'SampleZMprobaEllSub_visual.fits'

    master = Table.read(input)
    # default entries
    table = master['IAUNAME', 'PID', 'AID', 'RERUN',
                    'RUN', 'CAMCOL', 'FIELD', 'XPOS', 'YPOS'].copy()
    # change names for bdfitter
    table['IAUNAME'].name = 'NAME'
    table['PID'].name = 'PARENT_ID'
    table['AID'].name = 'ATLAS_ID'

    # rotation angle for SDSS PSF
    angle = []
    for row in master:
        name = 'data/sdss_field/frame-r-%06i-%i-%04i.fits' % (
            row['RUN'], row['CAMCOL'], row['FIELD'])
        angle.append(fits.getheader(name)['SPA'])
    table.add_column(Column(name='SPA', data=array(angle).astype(float)))

    # initial values
    Nrows = len(master)
    SER2_VAL = Column(
            name='SER2_VAL',
            data=array([[1, 10, 4, 0.7, 0, 0, 0, 0]*2]*Nrows).astype(float))
    SER2_FIX = Column(
            name='SER2_FIX',
            data=array([[0, 0, 0, 0, 1, 0, 0, 0]*2]*Nrows))
    table.add_columns([SER2_VAL, SER2_FIX])

    # sersic index
    table['SER2_VAL'][:,8+ip['n']] = 1.
    # position angle
    table['SER2_VAL'][:,ip['pa']] = deg2rad((master['SERSIC_PHI'].data + 90.) % 360.)
    table['SER2_VAL'][:,8+ip['pa']] = deg2rad((master['SERSIC_PHI'].data + 90.) % 360.)
    # axis ratio
    table['SER2_VAL'][:,ip['q']] = master['SERSIC_BA'].data
    table['SER2_VAL'][:,8+ip['q']] = master['SERSIC_BA'].data
    # effective radius
    table['SER2_VAL'][:,ip['Re']] = master['SERSIC_TH50'].data / 0.396
    table['SER2_VAL'][:,8+ip['Re']] = master['SERSIC_TH50'].data / 0.396 * 2.
    # center
    table['SER2_VAL'][:,ip['x']] = master['XCEN'].data
    table['SER2_VAL'][:,ip['y']] = master['YCEN'].data
    table['SER2_VAL'][:,8+ip['x']] = master['XCEN'].data
    table['SER2_VAL'][:,8+ip['y']] = master['YCEN'].data

    table.write(outname)






if __name__=='__main__':
    dvc_exp('sample.fits', 'fit/dvcexp/input_DVCEXP.fits')