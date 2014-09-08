""" script to prepare sdss field images for bdfitter

assumes following directory structure
data/
    nsa/
        pimages/ -- masking out stars
        images/ -- further deblending (child)
    sdss_field/ -- SDSS field images (eventually to be read from Peyton disk)
    sdss_ivar/ -- ivar for the entire field
    sdss_psf_meta/ -- SDSS psField files
    images/ -- final images are stored here
    ivar/ -- final ivar images are stored here
    hdr/ -- directory for montage hdr output
"""

import montage_wrapper as montage
from scipy.interpolate import interp2d
from astropy.table import Table
from astropy.io import fits
from astropy import wcs
from pylab import *
import logging, os

from mpltools import mkdir_p


formatter = logging.Formatter('%(asctime)s:%(levelname)s - %(message)s')
handler = logging.StreamHandler()
handler.setFormatter(formatter)
handler.setLevel(logging.INFO)

log = logging.getLogger(__name__)
log.addHandler(handler)
log.setLevel(logging.INFO)


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


def nearest_neighbor(x, y, z, xnew, ynew):
        xnew = atleast_1d(xnew)
        ynew = atleast_1d(ynew)
        if not len(xnew.shape)==1 and len(ynew.shape)==1:
            raise ValueError, "xnew and ynew should be 1-dim array"
        nx, ny = len(xnew), len(ynew)
        znew = zeros([ny, nx])
        for iy, ytemp in enumerate(ynew):
            for ix, xtemp in enumerate(xnew):
                iiy = abs(y-ytemp).argmin()
                iix = abs(x-xtemp).argmin()
                znew[iy, ix] = z[iiy, iix].copy()
        return znew


def nearest_neighbor_wcs(wcs, img, wcs_new):
    """
    Generate image according to new wcs
    no interpolation is done.
    """
    img_new = zeros([wcs_new.naxis2, wcs_new.naxis1])
    
    for irow, row in enumerate(img_new):
        for icol, pix in enumerate(row):
            ra, dec = wcs_new.wcs_pix2world(icol, irow, 1)
            iicol, iirow = wcs.wcs_world2pix(ra, dec, 1)
            iicol, iirow = int(iicol), int(iirow)
            try:
                img_new[irow, icol] = img[iirow, iicol].copy()
            except IndexError:
                img_new[irow, icol] = -1
    return img_new


def save_ivar(filter, run, camcol, field, datadir):
    """ Build and save inverse variance image """

    sdss_field_name = datadir + 'sdss_field/'+get_framename(run, camcol, field, filter)
    hdulist = fits.open(sdss_field_name)
    hdr, img = hdulist[0].header, hdulist[0].data
    nrowc = img.shape[0]

    # 2. read in sky, and interpolate to full image size; this returns a
    #    sky image the same size as the frame image, in units of counts
    sky = hdulist[2].data
    allsky = sky['allsky'][0]
    f = interp2d(arange(allsky.shape[1]), arange(allsky.shape[0]), allsky, kind='linear')
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

    ivar_out = sdss_field_name.replace('sdss_field', 'sdss_ivar')
    log.info("Save inverse variance image: %s", ivar_out)
    fits.PrimaryHDU(data=ivar, header=hdr).writeto(ivar_out, clobber=True)



master = Table.read('SampleZMprobaEllSub_visual.fits')

filter = 'r'
child_ext = 2
datadir = 'data/'
montage_silent = True
tsize = len(master)

# for igal, gal in enumerate(master):
if 1:
    gal = master[6]

    size = gal['SIZE']  # in deg
    ra, dec = gal['RA_1'], gal['DEC_1']
    run, camcol, field = gal['RUN'], gal['CAMCOL'], gal['FIELD']
    iauname, pid, aid = gal['IAUNAME'], gal['PID'], gal['AID']

    if not datadir.endswith('/'): datadir + '/'
    for subdir in ['sdss_ivar', 'images', 'ivar', 'hdr', 'diff', 'diffre', 'pimagere']:
        mkdir_p(datadir + subdir) 
    # input files
    field_name = datadir + 'sdss_field/'+get_framename(run, camcol, field, filter)
    ivar_name  = datadir + 'sdss_ivar/'+get_framename(run, camcol, field, filter)
    nsa_image_name = datadir + 'nsa/images/' + get_childname(iauname, pid, aid)
    nsa_parent_name = datadir + 'nsa/ivar/%s-parent-%s.fits.gz' % (iauname, pid)
    nsa_pimage_name = datadir + 'nsa/pimages/%s-pimage.fits.gz' % (iauname)
    # output files
    out_image_name = datadir + 'images/%s_cutout.fits' % (iauname)
    out_ivar_name = datadir + 'ivar/%s_ivar.fits' % (iauname)
    # tempfiles
    out_image_hdr =  datadir + 'hdr/%s_cutout.hdr' % (iauname)
    out_nsa_diff_name = datadir + 'diff/%s.fits' % (iauname)
    out_nsa_diff_re_name = datadir + 'diffre/%s.fits' % (iauname)
    out_nsa_pimage_re_name = datadir + 'pimagere/%s_nsa_pimage_re.fits' % (iauname)

    log.info("Process %s", iauname)
    if not os.path.isfile(ivar_name): save_ivar(filter, run, camcol, field)

    log.info("Cutout field image: %s", out_image_name)
    montage.mSubimage(field_name, out_image_name, ra, dec, size)
    log.info("Making header file: %s", out_image_hdr)
    hdr = montage.mGetHdr(out_image_name, out_image_hdr)

    # reproject everything to cutout field image
    log.info("Cutout ivar image: %s", out_ivar_name)
    montage.reproject(ivar_name, out_ivar_name, header=out_image_hdr, exact_size=True,
                        silent_cleanup=montage_silent)
    
    log.info("Save parent - child: %s", out_nsa_diff_name)
    hdu_child = fits.open(nsa_image_name)
    hdu_parent = fits.open(nsa_parent_name)
    fits.writeto(out_nsa_diff_name, data=hdu_parent[child_ext*2].data - hdu_child[child_ext].data,
                    header=hdu_child[2].header, clobber=True)        
    log.info("Reprojecting diff: %s", out_nsa_diff_re_name)
    montage.reproject(out_nsa_diff_name, out_nsa_diff_re_name, header=out_image_hdr, exact_size=True,
                        silent_cleanup=montage_silent)


    log.info("WCS reproject pimage: %s", out_nsa_pimage_re_name)
    img = fits.open(out_image_name, mode='update')
    wimg = wcs.WCS(img[0].header)
    imnaxis1, imnaxis2 = wimg.naxis1, wimg.naxis2
    imworldx = wimg.wcs_pix2world(arange(imnaxis1), zeros(imnaxis1), 1)[0]
    imworldy = wimg.wcs_pix2world(zeros(imnaxis2), arange(imnaxis2), 1)[1]

    pimg = fits.open(nsa_pimage_name)
    wpimg = wcs.WCS(pimg[0].header)
    pimnaxis1, pimnaxis2 = wpimg.naxis1, wpimg.naxis2
    pimworldx = wpimg.wcs_pix2world(arange(pimnaxis1), zeros(pimnaxis1), 1)[0]
    pimworldy = wpimg.wcs_pix2world(zeros(pimnaxis2), arange(pimnaxis2), 1)[1]
    
    pimg_new = nearest_neighbor_wcs(wpimg, pimg[0].data, wimg)
    log.info("written to %s", out_nsa_pimage_re_name)
    fits.writeto(out_nsa_pimage_re_name, pimg_new, img[0].header, clobber=True)

    # do pimage masking
    # find out target galaxy id
    pix_x, pix_y = wimg.wcs_world2pix(ra, dec, 1)  # pixel coord of target
    pix_x, pix_y = int(pix_x), int(pix_y)
    galid = pimg_new[pix_y, pix_x]
    log.info("Target galaxy id = %3i", galid)
    mask = ones([wimg.naxis2, wimg.naxis1])
    mask[where((pimg_new > 0) & ( abs(pimg_new-galid) > 0.5 ))] = 0.
    pmasked = img[0].data * mask

    # do deblending
    deblended = pmasked - fits.getdata(out_nsa_diff_re_name)
    img[0].data = deblended

    print 'update', img.filename()
    img.flush()
    img.close()

