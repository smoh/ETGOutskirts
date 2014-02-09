
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
- 'reject' sample: probaEll < 0.2 28,191
- 'select' sample: probaEll > 0.6 5407


For visual inspection, randomly select 250 galaxies from

- Sample 1 -- `rand_RS.cat`
- Sample 2 reject -- `rand_probaEll_reject.cat`
- Sample 2 select -- `rand_probaEll_select.cat`



* * *

We aim to constitute a representative(complete?) sample of local massive elliptical galaxies.
Our parent sample is the NASA-Sloan Atlas Catalog (citation) contatining galaxies
with $z < 0.055$ within the coverage of SDSS DR8 (cite). We select this catalog for
two main reasons:

    1) It complements bright nearby galaxies that may be in SDSS imaging footprint, but
    are eliminated due to difficulties in reduction (is this correct?), and adopts
    higher threshold for distriminating stars from galaxies.
    2) The catalog provides SDSS images with improved sky subtraction (cite), which is
    crucial in our analysis of faint outskirts of elliptical galaxies.




