#!/usr/bin/env python
# -*- coding: utf-8 -*-
###########################
#
# make_input.py
#
# this is an example of how to make input for fit_sample.pro
#
##############################
"""
Created on Thu Feb  6 11:39:36 2014

@author: clackner

program details
"""
from pylab import *
import sys, os, re
import argparse
import numpy as np
import pyfits

from astropy.table import Table, Column
import argparse
from configobj import ConfigObj

# parameter index dictionary
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

def single_DVC():
    """
    generate input table for single De Vauc profile

    inputtable : str
        input table name; should be part of NSA catalog in FITS format
    outname : str
        output table name
    """
    # this will be arguments at some point
    input = 'SampleZMprobaEllSub_visual.fits'
    outname = 'single_DVC/input_DVC.fits'

    master = Table.read(input)
    # default entries
    table = master['IAUNAME', 'PID', 'AID'].copy()
    # change names for bdfitter
    table['IAUNAME'].name = 'NAME'
    table['PID'].name = 'PARENT_ID'
    table['AID'].name = 'ATLAS_ID'

    Nrows = len(master)
    DVC_VAL = Column(
            name='DVC_VAL',
            data=array([[1, 10, 4, 0.7, 0, 0, 0, 0]]*Nrows).astype(float))
    DVC_FIX = Column(
            name='DVC_FIX',
            data=array([[0, 0, 1, 0, 1, 0, 0, 0]]*Nrows))
    table.add_columns([DVC_VAL, DVC_FIX])

    table['DVC_VAL'][:,pind['posang']] = deg2rad((master['SERSIC_PHI'].data + 90.) % 360.)
    table['DVC_VAL'][:,pind['ratio']] = master['SERSIC_BA'].data
    table['DVC_VAL'][:,pind['Reff']] = master['SERSIC_TH50'].data / 0.396
    table['DVC_VAL'][:,pind['center_x']] = master['XCEN'].data
    table['DVC_VAL'][:,pind['center_y']] = master['YCEN'].data

    table.write(outname)


def single_SER():
    """
    generate input table for single Sersic profile

    inputtable : str
        input table name; should be part of NSA catalog in FITS format
    outname : str
        output table name
    """
    # this will be arguments at some point
    input = 'SampleZMprobaEllSub_visual.fits'
    outname = 'single_SER/input_SER.fits'

    master = Table.read(input)
    # default entries
    table = master['IAUNAME', 'PID', 'AID'].copy()
    # change names for bdfitter
    table['IAUNAME'].name = 'NAME'
    table['PID'].name = 'PARENT_ID'
    table['AID'].name = 'ATLAS_ID'

    Nrows = len(master)
    SER_VAL = Column(
            name='SER_VAL',
            data=array([[1, 10, 4, 0.7, 0, 0, 0, 0]]*Nrows).astype(float))
    SER_FIX = Column(
            name='SER_FIX',
            data=array([[0, 0, 0, 0, 1, 0, 0, 0]]*Nrows))
    table.add_columns([SER_VAL, SER_FIX])

    table['SER_VAL'][:,pind['posang']] = deg2rad((master['SERSIC_PHI'].data + 90.) % 360.)
    table['SER_VAL'][:,pind['ratio']] = master['SERSIC_BA'].data
    table['SER_VAL'][:,pind['Reff']] = master['SERSIC_TH50'].data / 0.396
    table['SER_VAL'][:,pind['center_x']] = master['XCEN'].data
    table['SER_VAL'][:,pind['center_y']] = master['YCEN'].data

    table.write(outname)


def make_input_two():
    """
    Generate input table for de Vauc + Sersic
    """
    from astropy import cosmology
    cosmo = cosmology.FlatLambdaCDM(70., 0.3)
    master = Table.read('SampleZMprobaEllSub_visual.fits')

    # default entries
    table = master['IAUNAME', 'PID', 'AID'].copy()
    # change names for bdfitter
    table['IAUNAME'].name = 'NAME'
    table['PID'].name = 'PARENT_ID'
    table['AID'].name = 'ATLAS_ID'

    Nrows = len(master)
    DVC_VAL = Column(
            name='DVC_VAL',
            data=array([[1, 10, 4, 0.7, 0, 0, 0, 0]]*Nrows).astype(float))
    DVC_FIX = Column(
            name='DVC_FIX',
            data=array([[0, 1, 1, 0, 1, 0, 0, 0]]*Nrows))
    # SER_VAL = Column(
    #         name='SER_VAL',
    #         data=array([[1, 10, 4, 0.7, 0, 0, 0, 0]]*Nrows).astype(float))
    # SER_FIX = Column(
    #         name='SER_FIX',
    #         data=array([[0, 0, 1, 0, 1, 0, 0, 0]]*Nrows))
    table.add_columns([DVC_VAL, DVC_FIX])

    # posang
    table['DVC_VAL'][:, 7] = deg2rad((master['SERSIC_PHI'].data + 90.) % 360.)
    # centers
    table['DVC_VAL'][:,5] = master['XCEN'].data
    table['DVC_VAL'][:,6] = master['YCEN'].data
    # q_ba
    table['DVC_VAL'][:,3] = master['SERSIC_BA'].data

    
    # size using FP of Bernardi+03 and Petrosian mag
    petroflux = master['PETROFLUX'][:,4] * 3.631e-6   # Jy
    sigma = master['VDISP'] 
    L = petroflux * 4. * pi * (cosmo.luminosity_distance(master['Z_1']).value)  # erg/s/Hz
    log_R0_kpc = -2.34*log10(sigma) + 1.5*log10(L) + 15.592
    R0_pixels = 10**log_R0_kpc / (cosmo.angular_diameter_distance(master['Z_1']).value * 1000) \
                / 206265. / 0.396  # circularized effective radius in pixels
    table['DVC_VAL'][:,1] = R0_pixels

    table.write('input_DVC.fits')

    



def main():
    parser = argparse.ArgumentParser(description="Make input FITS table for bdfitter")    
    parser.add_argument('catalog', help="sample catalog file")
    parser.add_argument('config_file', help='configuration file')
    args = parser.parse_args()

    config = ConfigObj('test.cfg')

    # default entries
    master = Table.read(args.catalog)
    table = master['IAUNAME', 'PID', 'AID'].copy()
    # change names for bdfitter
    table['IAUNAME'].name = 'NAME'
    table['PID'].name = 'PARENT_ID'
    table['AID'].name = 'ATLAS_ID'

    # read profile info from config file
    for profile in config.sections:
        ncomponents = max(1, len(profile.sections))
        val = np.zeros(ncomponents*8)
        fix = np.zeros(ncomponents*8)

        # initial values from NSA 2-dim sersic fit result
        

        # fixed parameters
        if p.has_key('fix'):
            for parname, parvalue in p['fix']:
                val[pind[parname]] = parvalue
                fix[pind[parname]] = 1.
        # always fix shape to zero (*NOT USED*)
        val[pind['shape']] = 0.
        fix[pind['shape']] = 1.

    
    #now, if you don't want to use the default profiles with nothing held fixed
    #you need to set the input parameters here
    #the things to set are NAMEOFPROFILE_FIX, NAMEOFPROFILE_VAL
    #there should be arrays of the size of the output 
    #parameters (8x the number of profiles)
    #The name of the profile is whatever you want, but you need to give it
    #to fit_sample, so it can find the input parameters
    #XXX_FIX is a boolean array which =1 for fixed parameters
    #XXX_VAL is the input values, including fixed values
    #if you don't fix the initial position of flux, default values will 
    #be chosen by the code (based on image size)
    #the parameters (in order) are:
    #0--the surface brightness
    #1 --the half light radius (refers to the semimajor axis)
    #2 --Sersic index
    #3 --axis ratio (minor/major)
    #4 --shape of isophote, default is zero and fixed, 
    #   but if you specify XXX_FIX you need to fix the 4th parameter
    #5 --x coordinate of center of profile, be default multiple profiles 
    #    share a center
    #6 --y coordinate of center
    #7 --position angle (radians counterclockwise from x-axis)
    
    #first profile will be de Vauc with only the isophote shape fixed
    DVC_FIX = np.ndarray((len(names), 8))
    DVC_VAL = np.ndarray((len(names), 8))
    #second profile will be EXP+DVC with all the de Vauc. parameters except
    #the normalization fixed
    EXPDVC_FIX = np.ndarray((len(names), 16))
    EXPDVC_VAL = np.ndarray((len(names), 16))
    for i in range(len(names)):
        DVC_FIX[i,:] = [0,0,1,0,1,0,0,0]
        #you need to pick starting values for the size, and the axis ratio
        #starting values for the surface brightness and the position will
        #be set in the code, don't set the flux to zero, as the code
        #just rescales values (so can't rescale 0)
        DVC_VAL[i,:] = [1.0,10.,4.0,0.7,0.0,0.0,0.0,0.1]
        #for two component fits, always but the 'bulge' (smaller +higher sersic)
        #profile first
        EXPDVC_FIX[i,:] = [0,1,1,1,1,1,1,1,0,0,1,0,1,0,0,0]
        #you'll have to source the fixed values from somewhere, either
        #based on the radius, or based on previously fitting the profiles
        EXPDVC_VAL[i,:] = [1.0, 78.9, 4.0, 0.75, 0.0, 715.217, 725.011, 2.1801,
                            0.0, 100.0, 1.0, 0.75, 0.0, 0.0, 0.0, 2.2]
                            
         #it might be good to have different input scripts, for example
         #you could make 1 input file to fit the sersic/dvc profile and after
         #that runs, make another file for the two component fits, which 
         #takes the olds fits as inputs for the VALs arrays
    


    
    col4 = pyfits.Column(name='DVC_FIX', format='8D', array=DVC_FIX)
    col5 = pyfits.Column(name='DVC_VAL', format='8D', array=DVC_VAL)
    col6 = pyfits.Column(name='EXPDVC_FIX', format='16D', array=EXPDVC_FIX)
    col7 = pyfits.Column(name='EXPDVC_VAL', format='16D', array=EXPDVC_VAL)
    
    cols = pyfits.ColDefs([col1, col2, col3, col4,
                           col5, col6, col7])
    
    hdu = pyfits.new_table(cols)
    hdu.writeto('sample_input.fits', clobber=True)
    
   
    return 0




if __name__=='__main__':
    single_SER()
