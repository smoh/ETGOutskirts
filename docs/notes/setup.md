Claire's instructions on setting up her code
--------------------------------------------

In ~/.idlenv, put

```
#setup sloan stuff
source /u/dss/products/eups/bin/setups.sh
setup -v idlutils
setup -v photoop v1_9_11
setup -v idl v7_1
#I think version 8_2 also works, so feel free to change the above
#idl startup file
export IDL_STARTUP='/u/clackner/.idlstartup'

#add base idl library
export IDL_PATH=$IDL_DIR/lib:$IDL_PATH
export IDL_PATH=$IDL_PATH:$IDL_DIR/examples
export IDL_PATH=$IDL_PATH:+/u/clackner/pro
export IDL_DEVICE=X
```

In ~/.idlstartup
```
!except=2 ;show all exceptions where they occur
pref_set, 'idl_cpu_tpool_nthreads', 1, /commit ;do NOT use threading
set_plot, 'X'

;astrolib ;include Goddard libraries (alread in sdss)

pref_set, 'idl_rbuf_size', 1000, /commit ;increase buffersize to 1000
```

make_input.py which generates an example of the input file you'll need
and fit_sample.pro which is the actually fitter.


NSA Images
----------

Here's what I found on the NSA websites:
#NSA images
#all images are 5 images in bands: u,g,r,i,z,NUV(Galex),FUV(Galex)
#child image (with deblending)--this is what we want for fitting, I
refer to it as the image
http://sdss.physics.nyu.edu/mblanton/v0/detect/v0_1/RAh/pmXX/[IAUNAME]/atlases/PID/[IAUNAME]-[PID]-atlas-[CID].fits.gz

#parent images are are same dimension as child images. Images are odd
HDUs, Ivars are even HDUs, from this file, we just want the inverse
variance (
http://sdss.physics.nyu.edu/mblanton/v0/detect/v0_1/RAh/pmXX/[IAUNAME]/parents/[IAUNAME]-parent-[PID].fits.gz

#PSFs? Maybe, they look good, but I need to send an email to check
http://sdss.physics.nyu.edu/mblanton/v0/detect/v0_1/RAh/pmXX/[IAUNAME]/[IAUNAME]-[filter]-bpsf.fits.gz

If you are putting together a sample, the easiest way to organize things
is in 3 folders, one for the images (called images), one for the parent
images/inverse variance planes (called ivar) and one for the PSF (called
psf). I wrote input/output code that can get the right image for the
right band.

