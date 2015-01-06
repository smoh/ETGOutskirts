import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt
import numpy as np
from jug import TaskGenerator
import pybdfitter as bd

# Load my sample
s = bd.Sample('sample.fits', 'data', imgdir='admask')
s.add_model('deVExp', 'fit/deVExp/fit_deVExp.fits', 'fit/deVExp/models')
s.add_model('Ser', 'fit/Ser/fit_Ser.fits', 'fit/Ser/models')
# Set center to deVExp model
s.table['CENTER'] = np.vstack([s.models['deVExp'].p[0].x, s.models['deVExp'].p[0].y]).T

@TaskGenerator
def save_profile(ind):
    s.show_profile(ind, unit='mag', xunit='arcsec')
    plt.savefig("figures/profile/%03i_%s.png" % (ind, s.table['IAUNAME'][ind]))
    plt.close()

@TaskGenerator
def save_image(ind):
    s.show_image(ind)
    plt.savefig("figures/image/%03i_%s.png" % (ind, s.table['IAUNAME'][ind]), dpi=80)
    plt.close()

@TaskGenerator
def save_ellipse(ind):
    cx, cy = s.table['CENTER'][ind]
    outname = 'ellipse/%s.ell.dat' % (s.table['IAUNAME'][ind])
    try:
        bd.ellipse('data/admask/%s.fits' % (s.table['IAUNAME'][ind]),
                   outname,
                   x0=cx, y0=cy)
        fig = bd.EllipsePlot(outname)
        fig.setup_figure()
        plt.savefig(outname.replace('dat','png'))
        plt.close()
    except:
        print outname, 'failed'
    

for i in range(len(s.table)):
    save_profile(i)
    save_image(i)
    save_ellipse(i)

