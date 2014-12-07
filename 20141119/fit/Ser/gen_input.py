"""
Generate input file for single sersic model
"""

import numpy as np
from astropy.table import Table, Column
from astropy.io import fits

source_table = '../../testser/RAWFIT00000.00268.fits'
outname = 'input_ser.fits'

source = Table.read(source_table)

t_input = Table()
t_input.add_column(source['NAME'])

# put best-fit params of test Sersic fit as initial values
val = Column(
        name='SER_VAL',
        data=source['FIT_SER'].data)

fix = Column(
        name='SER_FIX',
        data=np.zeros((len(source), 8))
        )
fix[:, 4] = 1.  # do not fit shape of the isophotes
t_input.add_columns([val, fix])
t_input.write(outname)

