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


def save_ivar(filter, run, camcol, field):
    """ Build and save inverse variance image """

    sdss_field_name = 'data/sdss_field/'+get_framename(run, camcol, field, filter)
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
    print ivar_out, 'saved'
    hdu_ivar = fits.PrimaryHDU(data=ivar, header=hdr).writeto(ivar_out, clobber=True)



def reproject2child(filter, run, camcol, field, iauname, pid, aid, outname):
    """ reproject field image and ivar to match NSA child """
    
    sdss_field_name = 'data/sdss_field/'+get_framename(run, camcol, field, filter)
    sdss_ivar_name  = 'data/sdss_ivar/'+get_framename(run, camcol, field, filter)
    child = get_childname(iauname, pid, aid)
    child_ext = ['u','g','r','i','z'].index(filter)
    hdrout = 'data/hdr/' + child.replace('fits.gz','.hdr')
    montage.mGetHdr('data/nsa/images/'+child, hdrout, hdu=child_ext)  # extension important
    imgname = 'data/images/'+outname+'.fits'
    ivarname = 'data/ivar/'+outname+'.fits'
    montage.reproject(sdss_field_name, imgname, header=hdrout, exact_size=True)
    montage.reproject(sdss_ivar_name, ivarname, header=hdrout, exact_size=True)
    print imgname, 'saved'
    print ivarname, 'saved'


if 1:

    master = Table.read('SampleZMprobaEllSub_visual.fits')
    
    filter = 'r'
    child_ext = 2

for gal in master:   
    # gal = master[34]

    save_ivar(filter, gal['RUN'], gal['CAMCOL'], gal['FIELD'])
    reproject2child(filter, gal['RUN'], gal['CAMCOL'], gal['FIELD'],
                    gal['IAUNAME'], gal['PID'], gal['AID'], gal['IAUNAME'])


    imgname = 'data/images/' + gal['IAUNAME'] + '.fits'
    ra, dec = gal['RA_1'], gal['DEC_1']  # galaxy ra, dec in degrees
    img = fits.open(imgname, mode='update')
    wimg = wcs.WCS(img[0].header)
    imnaxis1, imnaxis2 = wimg.naxis1, wimg.naxis2
    imworldx = wimg.wcs_pix2world(arange(imnaxis1), zeros(imnaxis1), 1)[0]
    imworldy = wimg.wcs_pix2world(zeros(imnaxis2), arange(imnaxis2), 1)[1]

    pimage = 'data/nsa/pimages/'+'%s-pimage.fits.gz' % (gal['IAUNAME'])
    pimg = fits.open(pimage)
    wpimg = wcs.WCS(pimg[0].header)
    pimnaxis1, pimnaxis2 = pimg[0].header['NAXIS1'], pimg[0].header['NAXIS2']
    pimworldx = wpimg.wcs_pix2world(arange(pimnaxis1), zeros(pimnaxis1), 1)[0]
    pimworldy = wpimg.wcs_pix2world(zeros(pimnaxis2), arange(pimnaxis2), 1)[1]
    pimg_new = nearest_neighbor(pimworldx, pimworldy, pimg[0].data, imworldx, imworldy)
    pimg.close()

    child = fits.getdata('data/nsa/images/'+ \
        get_childname(gal['IAUNAME'], gal['PID'], gal['AID']), ext=child_ext)

    pix_x, pix_y = wimg.wcs_world2pix(ra, dec, 1)
    pix_x = int(pix_x)
    pix_y = int(pix_y)
    galid = pimg_new[pix_y, pix_x]

    pimg_new[where(pimg_new == galid)] = -1
    pimg_new[where(pimg_new > 0)] = 0
    pimg_new[where(pimg_new == -1)] = 1.
    pimg_new[where(child == 0)] = 0.

    img_masked = img[0].data * pimg_new
    img[0].data = img_masked.astype(float32)

    print 'update', img.filename()
    img.flush()
    img.close()
    pimg.close()








# if 0:
#     # use montage to reproject pimage
#     hdrout = imgname.replace('fits', 'hdr')
#     hdr = montage.mGetHdr(imgname, hdrout)
#     montage.reproject(pimage, 'pimage_reprojected.fits', header=hdrout, exact_size=True)



