
# The Fundamental Plane
## Bernardi et al. 2003 - I. The Sample

* Data
    - Surface photomerty by fitting de Vaucouleurs' profile accounting for the effects of seeing, atmospheric extinction, and Galactic extinction. Photon counts bined into radial + angular bins are fitted by model convolved with the PSF. From this,
        + minor-to-major axis ratio b/a
        + effective radius along major axis r_dev and total magnitude m_dev from deVauc
        + model magnitude from best-fit models in r* band
        + Petrosian magnitude (flux within 2 Petrosian radius)
        + _Petrosian radii?_
* Sample
    - About 9000 objects are selected based on concentration index, the ratio of likelihood of the de Vauc to the exponential model, and the PCA classification of spectra. Also requires that objects have z < 0.3, good spectra of S/N > 10. Galaxies with velocity dispersion less than 70 km/s are excluded.
    - resft-frame luminosities
        + cosmology
        + K-correction
    - rest-frame effective size
        + circular effective radius r_o = b/a * r_dev
        + _second correction I did not understand_
    - velocity dispersion
        + average of "Fourier-fitting" and "direct fitting" methods over 4000-7000 angstrom
        + error~0.02-0.06 dex depending on S/N
        + aperture correction
            * **assumes that ellipticals of different sizes have similar velocitys dispersion profile?**

## Bernardi et al. 2003 - III. FP            

* 9000 early-type galaxies in 0.01 < z < 0.3 (magnitude-limited)
* R_o $\propto$ sigma^1.49 I_o^-0.75 in r*
* little evolution
* Finding best-fit plane
    - "direct fit" and "orthogonal fit"
    - **maximum likelihodd approach** that accounts for evolution and selection effect v.s. least-square minimization
    - For FP as distance indicator, "direct fit" is more appropriate since it minimizes the scatter in size. For studies of stellar evolution, "orthogonal fit" makes more sense -- In their best-fit plane using maximum likelihood approach, the _intrinsic_ scatter in log r_o for direct and orthogonal fit is 0.088 and 0.094
    - 

## Saulder et al. 2013 - SDSS DR8

## La Barbera et al. 2010 - SPIDER
* SDSS DR6 5080 M_r < -25 ETGs with 0.05 < z < 0.095
* 

* * *

# Eigenthaler et al. 2013 - Age and metallicity gradients in fossil ellipticals
* Fossil galaxy groups: highly evolved system with central ellitical galaxy formed through multiple merging + group X-ray halo
* The radial gradient of _single-stellar population ages_ of six Fossil Central Galaxies is ~zero.
* The radial gradient of metallicities is negative with slope flatter than predicted by monolithic collapse, suggesting that these ellipticals are result of multiple major mergers.


* * *

SDSS imaging
- bands
    + u, g, r, i, z: 3560, 4680, 6180, 7500, 8870
- 0.396 arcsec/pixel
