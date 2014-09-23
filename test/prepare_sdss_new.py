""" script to prepare sdss field images for bdfitter

- Determine the desired cutout size as 1/sqrt(2) of pimage size so that cutout is
  covered regardless of relative rotation between cutout and pimage.
- If the nearest edge of the field image is smaller than the desired size, use
  this as cutout size.
- If the nearest edge of the field image is less than a quarter of the desired
  size, issue a flag.
- Prints index, IAU name, desired size (arcsec), actual size (arcsec), edge flag.
- Inverse variance images are first made from SDSS field images (with exactly same
  size), and chopped in the same dimension as cutouts.

- assumes the following directory structure
[datadir]/
    nsa/
        pimages/ -- labeled mask image
    sdss_field/ -- SDSS field images (eventually to be read from Peyton disk)
    sdss_psf_meta/ -- SDSS psField files
    sdss_ivar/ -- ivar for the entire field is stored here
    [outdir]/
        psf/
        cutout/ -- cutouts are stored here
        mask/ -- mask images are stored here
        fig/ -- cutout + mask image plots
        images/ -- final images are stored here
        ivar/ -- cutout ivar images are stored here
"""

import montage_wrapper as montage
from scipy import interpolate
from astropy.table import Table
from astropy.io import fits
from astropy import wcs
from pylab import *
import os

from mpltools import mkdir_p
from pybdfitter.sersic import bn, fn, nanomaggie2mag


def get_framename(run, camcol, field, filter):
    """
    Get sdss frame image name
    filter : one of 'u', 'g', 'r', 'i', 'z'
    """
    name = 'frame-%s-%06i-%i-%04i.fits' % (filter, run, camcol, field)
    return name

def get_childname(iauname, pid, aid):
    """
    Get NSA child image name
    """
    return '%s-%i-atlas-%i.fits.gz' % (iauname, pid, aid)


def sdss_gain(filter, camcol, run):
    """
    Return run-dependent gain for given filter and camcol

    filter : str, one of u, g, r, i, z
    camcol : int
    run : int

    Returns:
        gain : float
    """
    gainTable1 = Table.read("""
    camcol u g r i z
    1   1.62    3.32    4.71    5.165   4.745
    2   1.595   3.855   4.6 6.565   5.155
    3   1.59    3.845   4.72    4.86    4.885
    4   1.6 3.995   4.76    4.885   4.775
    5   1.47    4.05    4.725   4.64    3.48
    6   2.17    4.035   4.895   4.76    4.69
    """, format='ascii')
    gainTable2 = Table.read("""
    camcol u g r i z
    1   1.62    3.32    4.71    5.165   4.745
    2   1.825   3.855   4.6 6.565   5.155
    3   1.59    3.845   4.72    4.86    4.885
    4   1.6 3.995   4.76    4.885   4.775
    5   1.47    4.05    4.725   4.64    3.48
    6   2.17    4.035   4.895   4.76    4.69
    """, format='ascii')
    if run < 1100:
        return gainTable1[camcol-1][filter]
    else:
        return gainTable2[camcol-1][filter]

def sdss_dark_variance(filter, camcol, run):
    """
    Return run-dependent gain for given filter and camcol

    filter : str, one of u, g, r, i, z
    camcol : int
    run : int

    Returns:
        dark variance : float
    """
    dvTable1 = Table.read("""
    camcol  u   g   r   i   z
    1   9.61    15.6025 1.8225  7.84    0.81
    2   12.6025 1.44    1.00    5.76    1.0
    3   8.7025  1.3225  1.3225  4.6225  1.0
    4   12.6025 1.96    1.3225  6.25  9.61
    5   9.3025  1.1025  0.81    7.84    1.8225
    6   7.0225  1.8225  0.9025  5.0625  1.21
    """, format='ascii')
    dvTable2 = Table.read("""
    camcol  u   g   r   i   z
    1   9.61    15.6025 1.8225  7.84    0.81
    2   12.6025 1.44    1.00    6.25    1.0
    3   8.7025  1.3225  1.3225  4.6225  1.0
    4   12.6025 1.96    1.3225  7.5625  12.6025
    5   9.3025  1.1025  0.81    7.84    2.1025
    6   7.0225  1.8225  0.9025  5.0625  1.21
    """, format='ascii')
    if run < 1500:
        return dvTable1[camcol-1][filter]
    else:
        return dvTable2[camcol-1][filter]


def save_ivar(filter, run, camcol, field, datadir):
    """ Build and save inverse variance image """

    sdss_field_name = datadir+'/sdss_field/'+get_framename(run, camcol, field, filter)
    ivar_out = sdss_field_name.replace('sdss_field', 'sdss_ivar')  # output name
    
    # do nothing if the file exists
    if os.path.isfile(ivar_out):
        return ivar_out

    hdulist = fits.open(sdss_field_name)
    hdr, img = hdulist[0].header, hdulist[0].data
    nrowc = img.shape[0]

    # 2. read in sky, and interpolate to full image size; this returns a
    #    sky image the same size as the frame image, in units of counts
    sky = hdulist[2].data
    allsky = sky['allsky'][0]
    f = interpolate.interp2d(arange(allsky.shape[1]), arange(allsky.shape[0]), allsky, kind='linear')
    simg = f(sky['xinterp'][0], sky['yinterp'][0])  # sky image; same dimension as frame

    # 3. read in calibration, and expand to full image size; this returns
    #    a calibration image the same size as the frame image, in units of
    #    nanomaggies per count
    calib = hdulist[1].data
    cimg = calib * ones([nrowc, 1])

    gain = sdss_gain(filter, camcol, run)
    darkVariance = sdss_dark_variance(filter, camcol, run)
    dn = img/cimg + simg
    dn_err = sqrt(dn/gain + darkVariance)
    img_err = dn_err * cimg
    ivar = 1./img_err**2

    fits.PrimaryHDU(
        data=ivar, header=hdr).writeto(ivar_out, clobber=True, output_verify='ignore')
    return ivar_out


def save_psf(gal, filter, outdir, datadir='data'):
    """
    Save psf image using SDSS read_PSF
    
    """
    import subprocess
    filter_id = ['u','g','r','i','z'].index(filter) + 1
    psField_name = datadir+'/sdss_psf_meta/psField-%06i-%i-%04i.fit' % (
            gal['RUN'], gal['CAMCOL'], gal['FIELD'])
    psf_outname = outdir+'/%s.fits' % (gal['IAUNAME'])
    cmd = "read_PSF %s %i %f %f %s" % (
            psField_name, filter_id, gal['YPOS'], gal['XPOS'], psf_outname)
    # print cmd
    subprocess.call(cmd, shell=True)


def process_images(gal, filter, sizeopt=30., datadir='data', outdir='pipeline', dryrun=False):
    """
    Process galaxy images for bdfitter
    Prints the following info in a line:
        name
        desired cutout size in arcsec
        desired cutout size in pixels
        fraction actuall area / desired area of cutout

    TODO: complete docstring

    gal : dict-like
        must have these keys
        IAUNAME, RUN, CAMCOL, FIELD, ... 
    filter : str
        one of u, g, r, i, z
    sizeopt : str or float
        if string, it is the key for the desired cutout size in degrees
        if float, it is the cutoff surface brightness (mag/arcsec^2) using given Sersic model
        specified by these keys
            SERSIC_N : sersic index
            SERSIC_MU : surface brightness at Re
            SERSIC_TH50 : Re in arcsec
    """

    # Determine the size of the cutout image
    if type(sizeopt) is str:
        # size is explicitly specified in degrees
        cutout_size_deg = gal[sizeopt]
    elif type(sizeopt) is float:
        # use given Sersic best-fit to calculate the radius at which
        # the surface brightness (mag/arcsec^2) falls below *sizeopt*
        n = gal['SERSIC_N']  # sersic index
        mu_e = gal['SERSIC_MU']  # surface brightness at Re in mag/arcsec^2
        Re_pix = gal['SERSIC_TH50'] / 0.396  # Re in pixels
        x = (1. + log(10)/2.5/bn(n) * (sizeopt - mu_e))**n
        cutout_size_deg = 2*Re_pix * x * 0.396 / 3600.

    # Check if the cutout will large enough
    rectField = [0, 0, 2048, 1489]
    ce = r_[gal['XPOS'], gal['YPOS']]  # center pixel position in field image
    halfwidth = cutout_size_deg*3600/0.396/2  # in pixels
    rectCutout = [ce[0]-halfwidth, ce[1]-halfwidth, ce[0]+halfwidth, ce[1]+halfwidth]
    # overlapping fraction of desired cutout
    frac_overlap = (min(rectCutout[2], 2048) - max(rectCutout[0],0)) * \
            (min(rectCutout[3],1489) - max(rectCutout[1],0)) / (4*halfwidth**2)

    print gal['IAUNAME'], '%8.2f %8.2f %8.4f' % (cutout_size_deg*3600, 2*halfwidth, frac_overlap)

    if not dryrun:
        # Start processing
        # save psf
        save_psf(gal, filter, outdir+'/psf', datadir=datadir)
        # SDSS field image
        fn_field = datadir+'/sdss_field/'+get_framename(gal['RUN'], gal['CAMCOL'], gal['FIELD'], filter)
        # NSA mask image
        fn_pimage = datadir + '/nsa/pimages/'+'%s-pimage.fits.gz' % (gal['IAUNAME'])

        # cutout field image and open it
        fn_cutout = outdir+'/cutout/%s_cutout.fits'%(gal['IAUNAME'])
        montage.mSubimage(fn_field, fn_cutout, gal['RA_1'], gal['DEC_1'], cutout_size_deg)
        hdu_cutout = fits.open(fn_cutout)

        # make and cutout ivar image
        fn_field_ivar = save_ivar(filter, gal['RUN'], gal['CAMCOL'], gal['FIELD'], 'data/')
        fn_cutout_ivar = outdir+'/ivar/%s.fits' % (gal['IAUNAME'])
        montage.mSubimage(fn_field_ivar, fn_cutout_ivar, gal['RA_1'], gal['DEC_1'], cutout_size_deg)
        
        # resample pimage in cutout image coordinates using nearest neighbor method 
        # to preserve child id's exactly
        hdu_pimage = fits.open(fn_pimage)
        wcs_pimage = wcs.WCS(hdu_pimage[0].header, fix=False)
        x, y = meshgrid(arange(wcs_pimage.naxis1), arange(wcs_pimage.naxis2))
        wx, wy = wcs_pimage.all_pix2sky(x, y, 0)
        val = hdu_pimage[0].data

        wcs_cutout = wcs.WCS(hdu_cutout[0].header, fix=False)
        x1, y1 = meshgrid(arange(wcs_cutout.naxis1), arange(wcs_cutout.naxis2))
        wx1, wy1 = wcs_cutout.all_pix2sky(x1, y1, 0)
        pim = interpolate.griddata((wx.ravel(), wy.ravel()), val.ravel(), (wx1, wy1), method='nearest')
        # To prevent masks on the edge to extend unnecesarily due to nearest interpolation
        # reset values outside of pim to -1
        xout = (wx1<wx.min()) | (wx1>wx.max())
        yout = (wy1<wy.min()) | (wy1>wy.max())
        pim[xout | yout] = -1

        # find out child ID of the galaxy
        pix_x, pix_y = wcs_pimage.wcs_world2pix(gal['RA_1'], gal['DEC_1'], 0)#, ra_dec_order=True)
        gal_id = hdu_pimage[0].data[rint(pix_y), rint(pix_y)]

        # make a mask image with 0 for pixels that should be masked, 1 for others
        # to be multiplied to cutout image
        mask = where((pim != gal_id) & (pim != -1), 0., 1.)
        # write mask image
        hdu_mask = fits.HDUList([fits.PrimaryHDU(data=mask)])
        hdu_mask[0].header.update(wcs_cutout.to_header())  # put WCS header keywords
        hdu_mask.writeto(outdir+'/mask/%s_mask.fits' % (gal['IAUNAME']), clobber=True)
        
        # figure for visual check
        figure(figsize=(10,5))
        subplots_adjust(left=0.08, right=0.95)
        suptitle(' '.join([gal['IAUNAME'], fn_field.split('/')[-1]]), fontsize=10)
        subplot(121)
        imshow(hdu_cutout[0].data, norm=mpl.colors.LogNorm(), origin='lower', aspect='auto')
        subplot(122)
        # make a color map of fixed colors for mask image
        cmap = mpl.colors.ListedColormap(['black', 'white'])
        norm = mpl.colors.BoundaryNorm([0., .5, 1.], cmap.N)
        imshow(mask, interpolation='nearest', origin='lower', aspect='auto', cmap=cmap, norm=norm)
        savefig(outdir+'/fig/%s.png' % (gal['IAUNAME']), dpi=80)
        close()

        # save masked image
        hdu_cutout[0].data *= mask
        hdu_cutout.writeto(outdir+'/images/%s.fits' % (gal['IAUNAME']), clobber=True)


if __name__ == '__main__':
    nsa = Table.read('sample.fits')  # NSA table of the sample
    # add mu_e column
    Re = nsa['SERSIC_TH50'] / 0.396  # pixels
    ba = nsa['SERSIC_BA']
    Ie_avg = nsa['SERSICFLUX'][:,4] / (2.*pi*ba*(Re)**2)
    n = nsa['SERSIC_N']
    Ie = Ie_avg / fn(n)  # nanomaggie/pixel
    mu_e = nanomaggie2mag(Ie/.396**2)
    nsa['SERSIC_MU'] = mu_e

    outdir = 'pipeline'
    # make directory structure
    for subdir in ['psf', 'cutout', 'mask', 'fig', 'images', 'ivar']:
        mkdir_p(outdir+'/'+subdir)
    filter = 'r'
    for gal in nsa:
        process_images(gal, filter=filter, outdir=outdir)
    # process_images(nsa[0], filter='r', outdir=outdir)
