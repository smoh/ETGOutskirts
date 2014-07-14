# import matplotlib
# matplotlib.use('AGG')
from pylab import *
from pybdfitter.plot import Plotter, GalaxyImage, add_ellipse
from pybdfitter.sersic import NSersic, nanomaggie2mag, sersic2d
from pybdfitter import AsinhNorm, cosmo
from pybdfitter.fp import Bernardi03_r
import matplotlib.gridspec as gridspec
import time
import mpltools
import itertools
from astropy.table import Table
from jug import TaskGenerator

rc('font', family='monospace')
seterr(divide='warn')

rc('axes', grid=True)

# visual flag for deblended images
tvis = Table.read('visualflag.dat', format='ascii')
igood = where(tvis['vflag'] == '0')[0]


tin = Table.read('sample.fits')
t = Table.read('fit/dvc/RAWFIT00000.00268.fits')
tSER = Table.read('fit/ser/RAWFIT00000.00268.fits')
tDVCEXP = Table.read('fit/dvcexp/RAWFIT00000.00268.fits')

col_chisq_DVC = t['CHISQ_DVC'].data * t['DOF_DVC'].data
col_chisq_SER = tSER['CHISQ_SER'].data * tSER['DOF_SER'].data
col_chisq_DVCEXP = tDVCEXP['CHISQ_DVCEXP'].data * tDVCEXP['DOF_DVCEXP'].data

colsigma = tin['VDISP'].data
colz = tin['Z_1'].data
colkcorr = tin['KCORRECT'][:, 4].data
colext = tin['EXTINCTION'][:, 4].data
colDA = cosmo.angular_diameter_distance(colz).value*1e3/206265.

def R0_fp(mu_e_avg):
    I0 = (mu_e_avg - 10*log10(1+colz) - colkcorr - colext)/-2.5
    return 10**(Bernardi03_r.a*log10(colsigma) + Bernardi03_r.b*I0 + Bernardi03_r.c)

# <codecell>

colpDVC1 = NSersic(t['FIT_DVC'].data.T)
colR0_DVC_fp = R0_fp(colpDVC1[0].mu_e_avg)
colR0_DVC_fit = colpDVC1[0].Ro*0.396*colDA

colpSER = NSersic(tSER['FIT_SER'].data.T)
colR0_SER_fp = R0_fp(colpSER[0].mu_e_avg)
colR0_SER_fit = colpSER[0].Ro*0.396*colDA

colpDVCEXP = NSersic(tDVCEXP['FIT_DVCEXP'].data.T)
colR0_DVCEXP_DVC_fp = R0_fp(colpDVCEXP[0].mu_e_avg)
colR0_DVCEXP_DVC_fit = colpDVCEXP[0].Ro*0.396*colDA
colR0_DVCEXP_EXP_fp = R0_fp(colpDVCEXP[1].mu_e_avg)
colR0_DVCEXP_EXP_fit = colpDVCEXP[1].Ro*0.396*colDA

# scat_DVC = sqrt(mean((log10(colR0_DVC_fp)-log10(colR0_DVC_fit))**2))
# print scat_DVC
# scat_DVCEXP = sqrt(mean((log10(colR0_DVCEXP_DVC_fp)-log10(colR0_DVCEXP_DVC_fit))**2))
# print scat_DVCEXP

# @TaskGenerator
def do(i):
# for i in [22]: 
    name = t['NAME'][i]
    im = GalaxyImage('data/deblended/'+name+'.fits', NSersic(t['FIT_DVC'][i]),
                     ivar='data/ivar/'+name+'.fits')
    mDVC = GalaxyImage('fit/dvc/models/M'+name+'.fits', NSersic(t['FIT_DVC'][i]))
    mSER = GalaxyImage('fit/ser/models/M'+name+'.fits', NSersic(tSER['FIT_SER'][i]))
    mDVCEXP = GalaxyImage('fit/dvcexp/models/M'+name+'.fits', NSersic(tDVCEXP['FIT_DVCEXP'][i]))
    
    pDVC = NSersic(t['FIT_DVC'].data[i])
    pSER = NSersic(tSER['FIT_SER'].data[i])
    pDVCEXP = NSersic(tDVCEXP['FIT_DVCEXP'].data[i])
    
    chisq_DVC, chisq_SER, chisq_DVCEXP = col_chisq_DVC[i], col_chisq_SER[i], col_chisq_DVCEXP[i]

    t1 = time.time()
    figure(figsize=(11,8.5))
    
    # -------------------------------------------------------------------------
    # images
    # -------------------------------------------------------------------------
    gs = gridspec.GridSpec(
        1, 4,bottom=0.7, wspace=0.05, top=0.95, left=0.05, right=0.9)
    axim_list = [subplot(gs[0,j]) for j in arange(4)]
    cax = axes([0.92, 0.7, 0.02, 0.25])
    extent=[-im.xc*0.396,
            (im.data.shape[1]-im.xc)*0.396,
            -im.yc*0.396,
            (im.data.shape[0]-im.yc)*0.396]
    # set norm
    vmin = percentile(im.data.ravel(), 5)
    vmax = percentile(im.data.ravel(), 95)
    norm = AsinhNorm(vmin=vmin, vmax=vmax)
    cmap = 'gray'
    dic = {'data': (0, im, 0),
           'DVC': (1, mDVC, pDVC),
           'SER': (2, mSER, pSER),
           'DVCEXP': (3, mDVCEXP, pDVCEXP)}
    for k, (j, cim, cp) in dic.iteritems():
        if j == 0:
            aa = axim_list[j].imshow(
                cim.data, aspect='auto', origin='lower', extent=extent,
                cmap=cmap, norm=norm)#, norm=AsinhNorm())
            colorbar(aa, norm=norm, cax=cax)
        else:
            axim_list[j].imshow(
                im.data - cim.data, aspect='auto', origin='lower', extent=extent,
                cmap=cmap, norm=norm)
        axim_list[j].set_title(k)
        # effective radius ellipse
        if k is not 'data':
            for l in range(cp.nprofiles):
                add_ellipse(axim_list[j], 0., 0., #cp[l].x*.396, cp[l].y*.396,
                            cp[l].Re*.396, cp[l].Re*cp[l].q*.396, cp[l].pa, 
                            ec=['r','b'][l])
    for ax in axim_list[1:]: mpltools.hide_tick_labels(ax, 'y')

    # -------------------------------------------------------------------------
    # profiles
    # -------------------------------------------------------------------------
    gs = gridspec.GridSpec(
        4, 1, bottom=0.07, top=0.67, left=0.07, right=0.5, hspace=0.04)
    # surface brightness profiles --------------------------------------
    axp = subplot(gs[:2,:2])
    axp.errorbar(im.profile['r'], im.profile['sb'], im.profile['sberr'],
                 fmt='ko', capsize=0)
    axp.plot(mDVC.profile['r'], mDVC.profile['sb'],'r-', label='DVC')
    axp.plot(mSER.profile['r'], mSER.profile['sb'],'b-', label='SER')
    axp.plot(mDVCEXP.profile['r'], mDVCEXP.profile['sb'], 'm-', label='DVCEXP')
    # axp.axvline(pDVC[0].Re*.396, c='r', ls='--')
    # axp.axvline(pSER[0].Re*.396, c='b', ls='--')
    # axp.axvline(pDVCEXP[0].Re*.396, c='k', ls='--')
    # axp.axvline(pDVCEXP[1].Re*.396, c='k', ls='--')    
    # axp.set_xlabel("log (r / arcsec)")
    axp.set_ylabel("$\mu$ (mag/arcsec$^2$)")
    axp.set_ylim([min(30, max(im.profile['sb']*1.05)), min(im.profile['sb']*0.95)])
    # annotations
    s = [' DISK/TOT = %5.4f' % (tDVCEXP['FLUXRATIO_DVCEXP'][i]),
         '   DVC Re = %5.1f %5.2f' % (pDVC[0].Re, pDVC[0].Re*.396),
         '   SER Re = %5.1f %5.2f' % (pSER[0].Re, pSER[0].Re*.396),
         'DVCEXP Re = %5.1f %5.2f' % (pDVCEXP[0].Re, pDVCEXP[0].Re*.396),
         'DVCEXP Re = %5.1f %5.2f' % (pDVCEXP[1].Re, pDVCEXP[1].Re*.396)]
    axp.text(0.95, 0.95, '\n'.join(s),
             va='top', ha='right', transform=axp.transAxes)
    
    x, y = meshgrid(arange(im.data.shape[1]), arange(im.data.shape[0]))
    for k in range(pDVCEXP.nprofiles):
        m = GalaxyImage(sersic2d(x, y, pDVCEXP[k].p), pDVCEXP)
        axp.plot(m.profile['r'], m.profile['sb'], 'm:')
    mpltools.hide_tick_labels(axp, 'x')
    
    axp_res = subplot(gs[2,:2], sharex=axp)
    axp_res.fill_between(im.profile['r'], -im.profile['sberr'], im.profile['sberr'], color='0.5')
    axp_res.plot(im.profile['r'], im.profile['sb'] - mDVC.profile['sb'], 'ro-')
    axp_res.plot(im.profile['r'], im.profile['sb'] - mSER.profile['sb'], 'bo-')
    axp_res.plot(im.profile['r'], im.profile['sb'] - mDVCEXP.profile['sb'], 'mo-')
    axp_res.set_ylim(1., -1.)
    axp_res.set_yticks([-0.5, 0, 0.5])
    mpltools.hide_tick_labels(axp_res, 'x')
    # background level
    # median_nan = lambda x: median(x[~isnan(x)])
    # yy=median_nan(im.data[:100,:100].ravel()**2)
    # bg = nanomaggie2mag(sqrt(yy) /0.396**2)
    # axp.axhline(bg)
    
    # PA and axis ratio profiles --------------------------------------
    if 1:
        axpa = subplot(gs[3,:2], sharex=axp)
        axpa.set_xlabel('r (arcsec)')
        # axq = subplot(gs[5,:2])#, sharex=axp)
        axq = axpa.twinx()
        for cim, co in zip([im, mDVC, mSER, mDVCEXP], ['k', 'r', 'b', 'm']):
#             c = cim.calc_paq(90, 0.7)
            c = cim.paq
            axpa.errorbar(c['r'], c['pa'] % 180, c['er_pa'], c=co, capsize=0)#, marker='o')
            axq.errorbar(c['r'], c['q'], c['er_q'], c=co, capsize=0, ls='--')#, marker='o')
        axpa.set_ylim(0., 180.)
        axq.set_ylim(0.6, 1.1)
    
    t2 = time.time()
    print t2-t1

    # -------------------------------------------------------------------------
    # chi-square / FP
    # -------------------------------------------------------------------------
    gs = gridspec.GridSpec(2, 2, bottom=0.08, top=0.67, left=0.58, right=0.98)

    ax_fp = subplot(gs[:,:])
    ax_fp.plot(log10(colR0_DVC_fp[i]), log10(colR0_DVC_fit[i]),
               'ro', ms=15, alpha=0.5, label='DVC')
    ax_fp.plot(log10(colR0_SER_fp[i]), log10(colR0_SER_fit[i]),
               'bo', ms=15, alpha=0.5, label='SER')
    ax_fp.plot(log10(colR0_DVCEXP_DVC_fp[i]), log10(colR0_DVCEXP_DVC_fit[i]),
               'mo', ms=15, alpha=0.5, label='DVCEXP')
    ax_fp.plot(log10(colR0_DVCEXP_DVC_fp[i]), log10(colR0_DVCEXP_EXP_fit[i]),
               'o', c='m', alpha=0.5, ms=15)
    ax_fp.legend(loc='lower right', numpoints=1)
    ax_fp.set_xlabel('log(R0/kpc) FP')
    ax_fp.set_ylabel('log(R0/kpc) FITTING')

    xx = arange(-1., 1.5, 0.1)
    ax_fp.plot(xx, xx, 'k-', lw=1)

    dchi2_SER = col_chisq_SER[i] - col_chisq_DVC[i]
    dchi2_DVCEXP = col_chisq_DVCEXP[i] - col_chisq_DVC[i]
    s = [r'$\Delta \chi^2$    (SER-DVC) = %.2f' % (dchi2_SER),
         r'$\Delta \chi^2$ (DVCEXP-DVC) = %.2f' % (dchi2_DVCEXP)]
    ax_fp.text(
        -0.9, 1.4, '\n'.join(s),
        va='top', ha='left')

    # savefig('figs/%s.png' % (name))
    # close()



# for i in range(len(tin)):
#     do(i)
for i in [26]:
    do(i)
show()
