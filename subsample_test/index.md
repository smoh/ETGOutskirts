<style>
    a {color: #0033ff;}
    a:hover {background-color:#ffff66;}
    a:visited {color: #0033ff;}
</style>


## New subsample
 Bigger random selection from NSA
- [TOPCAT session](../sample.session.fits) from SampleZMprobaEll in NSA+MORPH catalog, assign a random number (0-1) and select random 10% with random number > 0.9 & VDISP>70 --> dubbed SampleZMprobaEll_sub N=359
- [testsample_ra_dec.txt](testsample_ra_dec.txt): visual classification using SDSS images
    * VELL column 1 if elliptical 0 if not
    * 136 excluded by visual (38%)
    * ra, dec files for SDSS ImgList: [testsample_ra_dec_v1.txt](testsample_ra_dec_v1.txt), [testsample_ra_dec_v0.txt](testsample_ra_dec_v0.txt)
    * testsample_vis.fits: choose only VELL==1

- `sample.fits`: These were combined with SampleZMprobaEll_visual. Excluding duplicates, total 269 galaxies
- download and prepare images for combined sample
- remove sdss_ivar to save space
* * *

# Subsample test

test sample of 49 galaxies

## with NSA psf
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
    - **masking**: use NSA mask
        + need the rotation angle to make north is up as in NSA or to reverse it --> achieved using montage
        + in child or parent? in child, pixels with flux=0 as in [this example](images/J141224_masking). [ivar](images/J141224_ivar.png) after [reprojecting](images/J141224_parent_child_cropped.png) using montage --> apparently background stars have also been removed and interpolated with sky noise..
        + use NSA pimages to mask out stars: [resulting images compared with NSA child](images/tile_compare_with_child.png), [J080720](images/J080720_compare_with_child.png), [J083629](images/J083629_compare_with_child.png)
2. need to generate inverse-variance images: follow [this](http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html) instruction by SDSS
3. PSF from SDSS at the given position in field image

### Deblending

- NSA images
    * mosaic
    * parent: cutout + masked (using pimage)
    * child: additional deblending    

## Preparing SDSS field images for bdfitter

* data directory structure
```
data/
    nsa/
        pimages/ -- masking out stars
        images/ -- further deblending
    sdss_field/ -- SDSS field images (eventually to be read from Peyton disk)
    sdss_psf_meta/ -- SDSS psField files
    images/
    ivar/
    psf/
```

## Fit results with SDSS images
* single de Vaucouleurs profile - [deblended](sdss_psf/dvc/deblended/summary.html)
* single Sersic profile
    - [deblended](sdss_psf/ser/deblended/summary.html), [with profile](sdss_psf/ser/deblended/summaryProfile.html)
    - [masked](sdss_psf/ser/masked/summary.html)
    - [raw](sdss_psf/ser/raw/summary.html)
* single de Vaucouleurs profile with Re fixed - 
* 2 Sersic profile - [deblended](sdss_psf/ser2/deblended/summary.html)
