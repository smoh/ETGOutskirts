;written to check consistency in results from this and python

;; 0. find filename of the frame file
framename = 'data/sdss_field/frame-r-005314-5-0233.fits'

;; 1. read in the FITS image from HDU0; the resulting image will be
;;    sky-subtracted as well as calibrated in nanomaggies/pixel
img= mrdfits(framename,0,hdr)
nrowc= (size(img,/dim))[1]

;; 2. read in sky, and interpolate to full image size; this returns a
;;    sky image the same size as the frame image, in units of counts
sky= mrdfits(framename,2)
simg= interpolate(sky.allsky, sky.xinterp, sky.yinterp, /grid)

;; 3. read in calibration, and expand to full image size; this returns
;;    a calibration image the same size as the frame image, in units of
;;    nanomaggies per count
calib= mrdfits(framename,1)
cimg= calib#replicate(1.,nrowc)

dn = img/cimg + simg
dn_err= sqrt(dn/4.725+0.81)
img_err= dn_err*cimg

; save all images
writefits, 'J141224_idl_sky.fits', simg
writefits, 'J141224_idl_calib.fits', cimg
writefits, 'J141224_idl_err.fits', img_err