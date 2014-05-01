""" script to prepare sdss field images for bdfitter

assumes following directory structure
[datadir]/
    nsa/
        pimages/ -- masking out stars
        images/ -- further deblending (child)
        iavr/ -- paret images
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
import uuid
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


def save_ivar(filter, run, camcol, field, datadir):
    """ Build and save inverse variance image """

    sdss_field_name = datadir+'/sdss_field/'+get_framename(run, camcol, field, filter)
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
    fits.PrimaryHDU(
        data=ivar, header=hdr).writeto(ivar_out, clobber=True, output_verify='ignore')



def reproject2child(filter, run, camcol, field, iauname, pid, aid, datadir):
    """ reproject field image and ivar to match NSA child """
    
    sdss_field_name = datadir+'/sdss_field/'+get_framename(run, camcol, field, filter)
    sdss_ivar_name  = datadir+'/sdss_ivar/'+get_framename(run, camcol, field, filter)
    child = get_childname(iauname, pid, aid)
    child_ext = ['u','g','r','i','z'].index(filter)
    hdrout = datadir+'/hdr/' + child.replace('fits.gz','.hdr')
    montage.mGetHdr(datadir+'/nsa/images/'+child, hdrout, hdu=child_ext)  # extension important
    imgname = datadir+'/raw/'+iauname+'.fits'
    ivarname = datadir+'/ivar/'+iauname+'.fits'
    montage.reproject(sdss_field_name, imgname, header=hdrout, exact_size=True,
                        silent_cleanup=True)
    montage.reproject(sdss_ivar_name, ivarname, header=hdrout, exact_size=True,
                        silent_cleanup=True)
    montage.mFixNan(ivarname, ivarname, nan_value=0.)
    print imgname, 'saved'
    print ivarname, 'saved'

    # put column position angle in the psf header
    # sdss_psf_name = datadir+'/sdss_psf/'+iauname + '-%s-psf.fits' % (filter)
    # angle = fits.getval(sdss_field_name, 'SPA', ext=0)
    # hdr = fits.getheader(sdss_psf_name)
    # hdr['SPA'] = angle



if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("catalog", type=str, help="NSA catalog FITS")
    parser.add_argument("datadir", type=str, help="data directory")
    parser.add_argument("-r", type=int, nargs=2,
                        help="-r low high; set to -1 for end")
    args = parser.parse_args()

    master = Table.read(args.catalog)    
    filter = 'r'
    child_ext = 2
    datadir = args.datadir

    # prepare directories
    for subdir in ['sdss_ivar', 'hdr', 'raw', 'masked', 'deblended', 'ivar']:
        mkdir_p(datadir + '/' + subdir)

    if args.r:
        low, high = args.r
        if high == -1:
            high = len(master)
    else:
        low, high = 0, len(master)

    for igal in range(low, high):
        run, camcol, field = master['RUN', 'CAMCOL', 'FIELD'][igal].data
        iauname, pid, aid = master['IAUNAME', 'PID', 'AID'][igal].data
        ra, dec = master['RA_1', 'DEC_1'][igal].data  # galaxy ra, dec in degrees
        print igal, iauname
        # save inverse variance of the frame image
        save_ivar(filter, run, camcol, field, datadir)
        # reproject frame/ivar to NSA child
        reproject2child(filter, run, camcol, field, iauname, pid, aid, datadir)

        # do masking and deblending of cutout frame image
        imgname = datadir + '/raw/' + iauname + '.fits'
        img = fits.open(imgname)
        wimg = wcs.WCS(img[0].header)
        imnaxis1, imnaxis2 = wimg.naxis1, wimg.naxis2
        imworldx = wimg.wcs_pix2world(arange(imnaxis1), zeros(imnaxis1), 1)[0]
        imworldy = wimg.wcs_pix2world(zeros(imnaxis2), arange(imnaxis2), 1)[1]

        # load mask image
        pimage = datadir + '/nsa/pimages/'+'%s-pimage.fits.gz' % (iauname)
        pimg = fits.open(pimage)
        wpimg = wcs.WCS(pimg[0].header)
        pimnaxis1, pimnaxis2 = pimg[0].header['NAXIS1'], pimg[0].header['NAXIS2']
        pimworldx = wpimg.wcs_pix2world(arange(pimnaxis1), zeros(pimnaxis1), 1)[0]
        pimworldy = wpimg.wcs_pix2world(zeros(pimnaxis2), arange(pimnaxis2), 1)[1]
        pimg_new = nearest_neighbor(pimworldx, pimworldy, pimg[0].data, imworldx, imworldy)
        pimg.close()
        pix_x, pix_y = wimg.wcs_world2pix(ra, dec, 1)
        pix_x = int(pix_x)
        pix_y = int(pix_y)
        galid = pimg_new[pix_y, pix_x]
        mask = zeros_like(pimg_new)  # mask to be multiplied
        mask[where(pimg_new == galid)] = 1.  # target
        mask[where(pimg_new == -1)] = 1.  # unmasked region

        # load extra sources to be deblended
        child = fits.getdata(datadir+'/nsa/images/'+ \
            get_childname(iauname, pid, aid), ext=child_ext)
        parent = fits.getdata(datadir+'/nsa/ivar/%s-parent-%s.fits.gz' % (iauname, pid),
                                ext=2*child_ext)
        extra = parent - child  # to be subtracted from image

        # save images
        img_masked = img[0].data * mask
        img[0].data = img_masked
        img.writeto(datadir+'/masked/'+iauname+'.fits', clobber=True,
                    output_verify='ignore') #, data=img_masked)
        img_deblended = (img[0].data - extra) * mask
        img_deblended[isnan(img_deblended)] = 0  # nan to zero
        img[0].data = img_deblended
        img.writeto(datadir+'/deblended/'+iauname+'.fits', clobber=True,
                    output_verify='ignore')

        img.close()
        pimg.close()


