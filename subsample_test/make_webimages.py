import jinja2
from glob import glob
from astropy.table import Table
import argparse

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
        <td><img width="80%" src="{{ '%s' % (plotdir) }}/{{ row['NAME'] }}.png"></td>
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

templatevars = {'table':result, 'profile':profiles[0], 'plotdir':args.plotdir}

with open(args.outname, 'w') as f:
    f.write(template.render(templatevars))