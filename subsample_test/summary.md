<style type="text/css">
    body {
        font-family: "HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif; 
        width: 800px;
        margin: 20px auto;
        text-align: justify;
        line-height: 1.5em;
    }
    table{
        border-collapse:collapse;
    }
    table, th, td {
        border: 1px solid black;
        padding: 5px;
    }
    a, a:visited {
        color: #6600FF;
    }
    a:hover {
        background-color: #99FFCC;
    }
</style>

2014-05-19

* [Notebook:TestChiSquare](http://nbviewer.ipython.org/url/www.astro.princeton.edu/~semyeong/etg/TestChiSquare.ipynb): F-test and AIC criterion applied


2014-05-08 Summary

## Sample

* I draw a bigger random selection from the parent sample being those
in the NSA catalog with the probability of begin Elliptical above 0.7 (Huertas-Company+11) and velocity dispersion > 70 km/s. Of 3861 galaxies in the parent sample, I randomly selected 10% (N=359) of them.
* To remove galaxies that are not early-type, I visually inspect SDSS color images, and excluded 38% (N=136). Combined with the smaller subsample I already had, I am left with 269 galaxies.
    - Table of RA, dec for those kept and excluded (may be used for SDSS SkyServer): [kept](testsample_ra_dec_v1.txt) [excluded](testsample_ra_dec_v0.txt)
* To select the cleanest images I could work with to e.g., test fitting algorithms and such, I visually [flagged](visualflag.dat) them for the following:
    - 0: OK (N=130)
    - 1: image truncated nearby or within the galaxy (N=21)
    - 2: nearby bright source (N=66)
    - 3: bad deblending due to astrometry (N=24)
    - 4: bad deblending (N=61)
* I was mainly concerned with coming up with the cleanest sample to test various things, so it is likely that I was too conservative on my definition of image being 'OK'.
* I checked that this subsample of cleanset images do not have a very different distribution in redshift or photometric parameter space (using those from NSA single Sersic fitting), although I tend to throw out more of larger/brighter galaxies [Notebook:ImageQuality](http://nbviewer.ipython.org/url/www.astro.princeton.edu/~semyeong/etg/ImageQuality.ipynb)

## Fitting

So far, I have tried the following models. All parameters (except indices) are free. For all models, the initial condition for PA, q, Re, x, y is taken from NSA single Sersic fitting unless indicated otherwise.

1. single De Vaucouleurs [images](fit/dvc/images.html)
    * initial conditions: fix Sersic index to 4
2. single Sersic [images](fit/ser/images.html)
    * initial conditions: start Sersic index at 4
3. De Vaucouleurs + Exponential [images](fit/dvcexp/images.html)
    * initial conditions
        + comp 1 : n=4(fixed), Re from NSA
        + comp 2 : n=1(fixed), Re twice of comp 1, Ie 0.1 of comp 1
4. De Vaucouleurs (FP) + Exponential
    * Initially I do a run of model 3 to get reasonable values for using the FP. From the second iteration, I use the sigma and surface brightness of the De Vaucouleurs component and the Bernardi+03 FP to calculate the new effective radius, and _fix_ it to that value. This is cell In[4] in [Notebook:FixReTest](http://nbviewer.ipython.org/url/www.astro.princeton.edu/~semyeong/etg/FixReTest.ipynb).

### Single component models

[Notebook:SingleComponent](http://nbviewer.ipython.org/url/www.astro.princeton.edu/~semyeong/etg/SingleComponent.ipynb)

* As a check, I plot the Fundamental Plane relation using the parameters from de Vaucouleurs, which are supposed to have been used to _define_ the FP in the first place. We do see the relation with little offset (0.0034dex) and the 1-sigma scatter of 0.089 dex (~23%) in logR0
* If I free the Sersic index,
    - the indices are typically higher than 4 for most galaxies
    - Re increases by a factor of 2-3

which is indicative that the model tries to put more light at the center and/or the outer part.

* Parameters from single Sersic model also roughly follow the FP relation, but with bigger offset (0.024 dex) and scatter (0.099 dex), and systematic trend in residual.

### Double component models

[Notebook:DVCEXP](http://nbviewer.ipython.org/url/www.astro.princeton.edu/~semyeong/etg/DVCEXP.ipynb)

* While the reduced chi^2 is not much of a discriminant, the chi^2 values are reduced by 10^3-10^4 in most of the cases, which, of course, is partly expected because of the increased number of free parameters.
* Compared to model 1, Re of dV component is smaller as the Exp component models the light from outer part.
* Interestingly, while dV component is generally smaller than that of model 1, this component seem to follow the same FP relation unlike the other, EXP component.
* Model 4 tries to use the FP to constrain the size of the inner dV component. I only tested this for one galaxy J123948.37. In [Notebook:FixReTest](http://nbviewer.ipython.org/url/www.astro.princeton.edu/~semyeong/etg/FixReTest.ipynb) I plot how the key parameters change for 20 iterations for dV component (red dot) and EXP component (blue line).
    - First, you can see that dV component virtually stays the same until ~10 iterations, while the EXP component wanders around.
    - The size of the inner component does not really change, but because I _fixed_ the Re completely, it does increase slightly (by sub-pixel size) during the first 10 loops, and shows runaway to some parameter space because of what the EXP does in the meanwhile. Thus, what I think we should do is to stop the iteration when Re of dV component falls within some range the expected value.
    - From the fact that the dV component of model 3 ("free" DVC+EXP) follows the FP relation, I expect that for many of the cases, the loop will not go on for long.
    - This particular galaxy is indicated in DVCEXP Notebook. You can see that the Re for model 1 and 3 are very similar, and the dV of model 3 follows the FP while the EXP component is ~0.5 dex larger than that expected from the FP.
    - The residual images of model 1 for this galaxy does seem to show some weak but discernable outer light above background noise while the inner part is under- subtracted, which is impoved in model 3. However, you can also see that the single Sersic model equally models the galaxy better.
    - I also note that for this galaxy, whether or not we would better fit the data with an extra outer component is not so clear in 1D profile.

## Models & Model Selection

### Models

* dV
* Sersic
* dV + exp
* dV + Sersic
* 2 Sersic
* 3 Sersic?
* Test for PSF effect
    - use mismatched PSF
    - mask out the central region

For each of the two component model, I could _fix_ the effective radius of the inner component so that it satisfies the FP with some allowed scatter with similar iterations like model 4.

### Model Selection Criteria

* residual image
* 1D profiles of surface brightness, ellipticity, PA
* Quantitative measure of goodness-of-fit
    -  Bayesian Information Criterion (Schwarz 1978): assumes the data points are independent and identically distributed.
    -  Aiaike Information Criterion
    -  F-test

