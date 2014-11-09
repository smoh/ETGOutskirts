import os
import numpy as np
from astropy.io import fits
from astropy.table import Table

table = Table.read('../s.fits')
Nsample = 100

ff = open('runall.cl', 'w')
for name, center in table['IAUNAME', 'CENTER']:
    print name
    imgname = '../pipeline/images_admask/%s.fits' % (name)
    ivarname = '../pipeline/ivar/%s.fits' % (name)
    outdir = name[:5]
    os.mkdir(outdir)

    image = fits.getdata(imgname)
    ivar = fits.getdata(ivarname)
    badpixels = image == 0

    # write a cl scri
    clscript = """
cd /u/semyeong/projects/ETGOutskirts/test2/ellipse_resample
stsdas
isophote

unlearn ellipse geompar controlpar samplepar


# geompar
geompar.x0 = {2[0]:f}
geompar.y0 = {2[1]:f}
geompar.step = 2
geompar.linear = yes

ellipse.verbose = no

tprint.pwidth=INDEF

# run
for (i = 0; i < {1:d}; i +=1) {{
    ellipse("{0:s}/"//i, "{0:s}/"//i//".out")
    tprint("{0:s}/"//i//".out", >"{0:s}/"//i//".dat")
}}
""".format(outdir, Nsample, center)
    with open("run%s.cl" % (outdir), 'w') as f:
        f.write(clscript) 
    ff.write("cl < run%s.cl\n" % (outdir))

    # resample the image
    for i in range(Nsample):
        newimage = image + np.random.normal(0, 1./np.sqrt(ivar))
        newimage[badpixels] = 0
        hdu = fits.PrimaryHDU(newimage)
        fits.HDUList([hdu]).writeto(outdir+'/%i.fits' % (i))

ff.close()
