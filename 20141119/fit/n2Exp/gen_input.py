"""
Generate input file for n=2 + Exp model
"""

import numpy as np
from astropy.table import Table, Column
from astropy.io import fits

source_table = '../deV/fit_deV.fits'
outname = 'input_n2Exp.fits'

source = Table.read(source_table)

t_input = Table()
t_input.add_column(source['NAME'])

val = Column(
        name='n2Exp_VAL',
        data=np.zeros([len(source), 16]))
# put best-fit params of single deV fit as initial values for the deV comp
val[:, :8] = source['FIT_DEV']
val[:, 2] = 2.  # n_inner = 2

val[:, 8+0] = val[:, 0]*0.5
val[:, 8+1] = 2*val[:, 1]
val[:, 8+2] = 1.  # sersic index
val[:, 8+3] = 0.8  # axis ratio
val[:, 13:15] = val[:, 5:7]  # center

fix = Column(
        name='n2Exp_FIX',
        data=np.zeros((len(source), 16)))
# inner comp ---------------------------
fix[:, 2] = 1.  # sersic index
fix[:, 4] = 1.  # do not fit shape of the isophotes
# outer comp ---------------------------
fix[:, 8+2] = 1.  # sersic index
fix[:, 8+4] = 1. # do not fit shape of the isophotes

t_input.add_columns([val, fix])
t_input.write(outname)

