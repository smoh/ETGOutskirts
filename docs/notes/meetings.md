
## 2014-02-10 Skype w/ Claire on setting things up

- ivar: inside parent images
- take a closer look at `make_input.py (I could write my own version of this)`, 
  `fit_sample.pro`, `sersic.pro`
- she will be working on giving upper/lower limit to parameters (rather than just fixed). This is not high priority at the moment.
- should take a look at condor system at Peyton
- the main dependency of the whole code is idlutils.. but maybe others, too.
- Let's start with a small sample of 10-50 that covers a range of redshift, sigma, $M_r$, $R_e$

## 2014-03-06 Jenny

In any case, we need to infer size from dispersion. Try two things:

1. sigma --(Faber-Jackson relation)--> L --(Kormendy relation)--> size
2. define our own size-velocity dispersion relation using SDSS deVac model mag (see van der Wel et al. 2008 -- large scatter): median / some smaller size constraint for compact objects

## 2014-03-13 Claire

Discussed three issues in using the Fundamnetal Plane

1. the FP is defined with cirularized effective radius: fix it by initial axis ratio -- the difference between this and actually fixing the circularized radius(a combination of effective radius in semi-major axis, and axis ratio) should not be significant with b/a > 0.6
2. velocity dispersion aperture correction itself has effective radius in it -- negligible
3. K-correction: we can apply K-correction given in NSA catalog for all pixels

Point-Spread Function issue

- PSF images given in NSA are sometimes largely asymmetric, irregular, and messy -- possibly just some weighted average of adjacent stars..
- switch to SDSS PSF: given run, camcol, field, x, y position within the field, use SDSS data(in the disk at Peyton) and Robert's routine to calculate PSF

To do:

1. Check that not many of the galaxies in my test sample falls on field edges. --> quite a few does spread over two frames, but not many are cut down right in the middle.
2. Check that those run, camcol, field.. in parent/child image headers matches those values in NSA catalog --> yes, of the multiple images in header that went into mosaick, the one that best containes the galaxy is in the catalog.
3. Check if NSA images are re-oriented from SDSS images. --> yes, as such that north it up.


## 2014-03-17 Jenny and Claire

If SDSS DR10 has images with sky subtracted in the same way as NSA (Blanton et al. 2011), there is no reason to use mosaicked NSA images as this causes the issue of how to deal with PSF of mosaicked images.

We will use:

1. SDSS DR 10 sky-subtracted frames
    - cutout: use XPOS, YPOS in NSA catalog, and cutout square of, say, 5 SDSS effective radius around the galaxy
    - masking: use NSA mask
        + need the rotation angle to make north is up as in NSA or to reverse it --> Robert
        + in child or parent? check
2. need to generate inverse-variance images: Claire has a code for this. may need to download a metadata file for this.
3. PSF from SDSS at the given position in field image

TODO:

1. Confirm with MB that DR10 images used sky-subtraction method of Blanton et al. 2011 --> yes, correct
2. Ask RL about how to rotate the field images such that north is up to math NSA masks.
3. Claire will send the code to generate ivar images.


