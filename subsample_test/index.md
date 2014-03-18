Subsample test
==============

test sample of 49 galaxies

with NSA psf
------------
- [single DeVauc](single_DVC/summary.html)
- [single Sersic](single_SER/summary.html)
- [NSA mosaicked images](images/subsample_image.png)
- NSA PSF images in [linear](images/subsample_psf_linearscale.png), [logscale](images/subsample_psf_logscale.png), and [zscale](images/subsample_psf_zscale.png)

- NSA Mosaicking
    - Even though a galaxy seems to be well-contained within a single field, two or more adjacent fields may have been used to meet certain size criteria of mosaic around a galaxy.
    - While quite a few spread over two frames, not many are cut down right in the middle.
    - These images are rotated from SDSS field images so that north is up.
    - In header, the field images that went into the mosaick are recorded as `idR-{RUN#}-{FILTER}{COLUMN}-{FIELD}.fit` with header cards `FILE001, FILE002, ...`
- NSA PSFs are just average of stars in mosaicked image (i.e., could be from multiple field images)

- [SDSS psf](images/subsample_sdss_psf_zscale.png)

## Using SDSS field images

1. SDSS DR 10 sky-subtracted frames
    - cutout: use XPOS, YPOS in NSA catalog, and cutout square of, say, 5 SDSS effective radius around the galaxy
    - masking: use NSA mask
        + need the rotation angle to make north is up as in NSA or to reverse it --> Robert
        + in child or parent? in child, pixels with flux=0 as in [this example](images/J141224_masking). [ivar](images/J141224_ivar.png) does not contain information about mask.)
2. need to generate inverse-variance images
3. PSF from SDSS at the given position in field image

