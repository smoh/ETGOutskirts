import jinja2
from glob import glob
from astropy.table import Table, Column
import argparse
from numpy import sign

parser = argparse.ArgumentParser(
            description="make html page of fitting result")
parser.add_argument("result", type=str, help="fit result (RAWXXX.fits)")
parser.add_argument("plotdir", type=str, help="directory of png images")
parser.add_argument("outname", type=str, help="output html filename")
args = parser.parse_args()


# open the template from file
# templateLoader = jinja2.FileSystemLoader(searchpath='./')
# templateEnv = jinja2.Environment(loader=templateLoader)
# template = templateEnv.get_template('webimages.jinja')

template = jinja2.Template(
"""
<!DOCTYPE html>
<html>
<head>
</head>
<body>
<table border="1">
    {% for row in table %}
    <tr>
        <td width="30%">
            <table>
                
                <tr>
                    <td>NAME:</td>
                    <td>{{ row['NAME'] }}</td>
                </tr>
                <tr>
                    <td>{{'CHISQ_%s' % (profile)}}:</td>
                    <td>{{ row['CHISQ_%s' % (profile)] }}</td>
                </tr>
                <tr>
                    <td>{{'STAT_%s' % (profile)}}:</td>
                    <td>{{ row['STAT_%s' % (profile)] }}</td>
                </tr>
                <tr>
                    <td>REFF{{ '*' if row['%s_FIX' % (profile)][1] else '' }}:</td>
                    <td>{{ row['FIT_%s' % (profile)][1] }}</td>
                </tr>
                <tr>
                    <td>SERSIC_N{{ '*' if row['%s_FIX' % (profile)][2] else '' }}:</td>
                    <td>{{ row['FIT_%s' % (profile)][2] }}</td>
                </tr>
                <tr>
                    <td>Q_BA:</td>
                    <td>{{ row['FIT_%s' % (profile)][3] }}</td>
                </tr>
                
            </table>
        </td>
        <td>
            <img width="80%" src="{{ '%s' % (plotdir) }}/{{ row['NAME'] }}.png">
            <img src="http://skyservice.pha.jhu.edu/DR9/ImgCutout/getjpeg.aspx?ra={{radec[loop.index0][0]}}&dec={{radec[loop.index0][1]}}">
        </td>
    </tr>
    {% endfor %}
</table>
</body>
</html>
""")

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


templatevars = {'table':result, 'profile':profiles[0], 'plotdir':args.plotdir,
                'radec':radec}

with open(args.outname, 'w') as f:
    f.write(template.render(templatevars))
