import jinja2
from glob import glob
from astropy.table import Table, Column
import argparse
import os
from numpy import sign

template = jinja2.Template("images.html")

parser = argparse.ArgumentParser(
            description="make html page of fitting result")
parser.add_argument("result", type=str, help="fit result (RAWXXX.fits)")
parser.add_argument("plotdir", type=str, help="directory of png images")
parser.add_argument("outname", type=str, help="output html filename")
parser.add_argument("--opt", type=str, help="SDSS navigate options", default='G')
parser.add_argument("--sdss", action='store_true', help='get sdss cutout image?')
args = parser.parse_args()

plotdir_rel = os.path.relpath(args.plotdir, os.path.dirname(args.outname))

result = Table.read(args.result)
# find out the profile name
profiles = [x.split('FIT_')[-1] for x in result.colnames if 'FIT_' in x]
if len(profiles) > 1:
    raise 'More than one profile found'
# ra, dec from IAUNAME
def iauname2radec(name):
    sig = '-' if '-' in name else '+'
    rastr = name.split(sig)[0][1:]
    decstr = name.split(sig)[-1]
    rahr = float(rastr[:2])
    ramin = float(rastr[2:4])
    rasec = float(rastr[4:])
    radeg = rahr*15. + ramin*15./60. + rasec*15./3600.
    
    decde = float(decstr[:2])
    decmin = float(decstr[2:4])
    decsec = float(decstr[4:])
    decdeg = decde + decmin/60. + decsec/3600.
    decdeg = decdeg*sign(float(sig+'1'))

    return radeg, decdeg
radec = [iauname2radec(name) for name in result['NAME']]


templatevars = {'table':result, 'profile':profiles[0], 'plotdir':plotdir_rel,
                'radec':radec, 'opt':args.opt, 'sdss':args.sdss}

with open(args.outname, 'w') as f:
    f.write(template.render(templatevars))
