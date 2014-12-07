"""
Generate input file for single sersic model
"""

import numpy as np
from astropy.table import Table, Column
from astropy.io import fits

source_table = '../sample.fits'
outname = 'input_ser.fits'


source = Table.read(source_table)
# get center from cutout log
tlog = Table.read('../cutout.log', format='ascii')

t_input = Table()
t_input.add_column(source['IAUNAME'])
t_input['IAUNAME'].name = 'NAME'

val = Column(
        name='SER_VAL',
        data=np.zeros((len(source), 8))
        )
val[:, 0] = .1  # Ie, only needs to be non-zero
val[:, 1] = 10.  # half-light radius
val[:, 2] = 4.  # sersic index
val[:, 3] = 0.8  # minor/major axis ratio
val[:, 5] = tlog['col5'].data  # center x
val[:, 6] = tlog['col6'].data  # center y
val[:, 7] = 0.  # position angle (radians, ccw from x-axis)

fix = Column(
        name='SER_FIX',
        data=np.zeros((len(source), 8))
        )
fix[:, 4] = 1.  # do not fit shape of the isophotes
t_input.add_columns([val, fix])
t_input.write(outname)

