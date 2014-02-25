import jinja2
from glob import glob
from astropy.table import Table

templateLoader = jinja2.FileSystemLoader(searchpath='./')
templateEnv = jinja2.Environment(loader=templateLoader)
template = templateEnv.get_template('webimages.html')


images = glob('./out/*.png')

result = Table.read('./out/RAWFIT00000.00004.fits')


templatevars = {'title':'Test', 'table':result}

with open('test.html', 'w') as f:
    f.write(template.render(templatevars))