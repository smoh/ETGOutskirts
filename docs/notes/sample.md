
## Zhu, Blanton & Moustakas 2010
Sample

- parent sample: NYU-VAGC(all galaxies in SDSS with z<0.05) DR6 version 77,149 + RC3(galaxeis in SDSS imaging but without spectroscopy) 10,474
- red sequence: luminosity-dependent color cut (includes blue ellipticals) 37,026
    * 32,726 with SDSS spec
- preselect spheroidals:
    * "bulge + disk" model fitting of 1D surface brightness profile
    * B/(B+D) > 0.7 and ellipticity < 0.6 (equivalent to b/a > 0.4)
    * "featurelessness": ?
    * velocity dispersion > 70 km/s considering SDSS instrumental dispersion (22,621/32,726)
    * left with 2648
- visual inspection to rule out contaminants including
    * bulge-dominated SB0
    * galaxies with faint dust lanes
    * S0 with weak spiral structure
- final sample: 1923 with velocity dispersion > 70 km/s + 430 without SDSS spec


## Bernardi+13

- Bayesian Automated morhpological classifier (Huertas-Company+11)
- Nair & Abraham 2010 provides visual morph classification 0.01 < z < 0.1


## Experiment following Zhu+10

2014-02-07

Parent sample: NSA catalog 145,155

**Sample RS**

Construct a sample following similar procedures of Zhu et al. 2010

- Z < 0.05 123,396
- velocity dispersion > 70 km/s 42,101
- red sequence: $M_g - M_i > -0.05 \times (M_r + 16.0) + 0.65$ 
- ellipticity: $\epsilon = 1-b/a$
    + b/a from Stokes parameters at some radius
    + b/a from 2D Sersic fit <-- used
---> 22,378
NOTE: By using 'VDISP > 70', I restricted the sample to those with VDISP measurement.

Matched with Huertas-Company+11 (sep < 1 arcsec)

* Of 21,587 matched with HC+11, 13300 have probability of being elliptical less than
0.2 (~60%), i.e., this sample is obviously highly contaminated. The most stringent
constraint on the sample actually comes from surface brightness fitting.
Only 5007 have probability of being elliptical more than 0.6. These two subgroups
are contained within 'reject', and 'select' samples of Sample 2 -- all the conditions
above does not add anything to probaEll condition.

**Sample probaEll**

Use Bayesian Automated Classification given by Huertas-Company et al. 2011

- match with HC+11 107,125
- z < 0.05 88,912
- VD > 70 km/s 39528
- 'reject' sample: probaEll < 0.2 30829
- 'select' sample: probaEll > 0.6 5104

NOTE: By using 'VDISP > 70', I restricted the sample to those with VDISP measurement.
NOTE: The catalog of Huertas-Company+11 is based on SDSS DR7 spectroscopy catalog

For visual inspection, randomly select 250 galaxies from

- Sample 1 -- `rand_RS.cat`
    * clearly contains too many non-elliptical galaxies (>50%)
- Sample 2 reject -- `rand_probaEll_reject.cat`
    * ellipticals: 0 0 0 0 0 0 0 1 1 0 (in chunks of 25 galaxies using SDSS skyserver image list)
    * Practically no ellipticals are included.
- Sample 2 select -- `rand_probaEll_select.cat`
    * non-ellipticals: 4 5 4 4 5 11 6 2 3 6 --> 20 % contamination
    * contaminants: rings, bars, weak spiral structures, blue nuclear

NOTE on catalog:
    - VDISP values are not always from SDSS even in case the value exists in SDSS source catalog
    - The master catalog does not contain error on VDISP


* * * 

2013-02-14

Use Bayesian Automated Classification given by Huertas-Company et al. 2011

NSA All: 145,155
SampleZ: z < 0.05 123,396
hasMorph: matched with HC11 123,396 (sep < 1 arcsec, 90% within 0.4 arcsec)
SampleZM: z < 0.05 & has morph 88,912
SampleZMprobaEll: probaEll > 0.7 3861

* The entire SampleZM shows bimodal distribution of probability of being E/S0 or Elliptical,
with most of them below 0.2, and another peak around ~0.8-0.9. Note also that in testing
their classification with visually classified sample of Nair & Abraham 2010, they find
those visually classified as ellipticals have probaEll distribution peaked at ~0.8.
* One key concern in using Huertas+Company morphology catalog is since they base their
catalog on SDSS DR7 spectroscopy catalog, we may be biased against low-redshift, bright
galaxies(upper left corner of redshift versus $M_r$ plot).
* Elliptical galaxies selected with probaEll > 0.7 threshold ranges velocity dispersion
from 50-300 km/s, $M_r$ from -17 to -23, and effective radius from 1 to 20 arcsec (1 to 20 kpc).


SampleZMprobaEllSub: ($0 % 100 == 1 & _7 & Z_1>0.02) | ($0 % 4 == 0 & _7 & Z_1<0.02)
(7 = SampleZMprobaEll)

Subsample of N=75 with reasonable ranges in $M_r$, R_eff, and VDISP saved as SampleZMprobaEllsub.fits


* * *

We aim to constitute a representative(complete?) sample of local massive elliptical galaxies.
Our parent sample is the NASA-Sloan Atlas Catalog (citation) contatining galaxies
with $z < 0.055$ within the coverage of SDSS DR8 (cite). We select this catalog for
two main reasons:

    1) It complements bright nearby galaxies that may be in SDSS imaging footprint, but
    are eliminated due to difficulties in reduction (is this correct?), and adopts
    higher threshold for discriminating stars from galaxies.
    2) The catalog provides SDSS images with improved sky subtraction (cite), which is
    crucial in our analysis of faint outskirts of elliptical galaxies.




