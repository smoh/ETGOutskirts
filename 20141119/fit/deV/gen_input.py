"""
Generate input file for single deV model
"""

import numpy as np
from astropy.table import Table, Column
from astropy.io import fits

source_table = '../Ser/fit_Ser.fits'
outname = 'input_deV.fits'

source = Table.read(source_table)

t_input = Table()
t_input.add_column(source['NAME'])

# put best-fit params of single Sersic fit as initial values
val = Column(
        name='deV_VAL',
        data=source['FIT_SER'].data)
val[:, 2] = 4.  # sersic index

fix = Column(
        name='deV_FIX',
        data=np.zeros((len(source), 8)))
fix[:, 2] = 1.  # sersic index
fix[:, 4] = 1.  # do not fit shape of the isophotes
t_input.add_columns([val, fix])
t_input.write(outname)

