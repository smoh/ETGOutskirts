2014-09-05 Fri meeting with Jenny

1. Check FP, color of two component galaxies
 - color is not yet available (only r-band catalog present)
  - and check these galaxies with claire's fitter

2. the general trend with two-component fitting falls into three categories
  - *smooth* single Sersic profile fits well
    - Exp component never dominates, and is basically non-existent
  - when SerExp model equally fits the data well (or better), the 'bulge' sersic index decreases, and the component dominance switches at some 'deflection radius'
  - the radius of Exp is so large that it maybe sky???

3. sample
  - there are tens of thousands of galaxies without H11 morphology (~ without SDSS spec) -- do we exclude these initially, and have another ``special'' sample?
  - clear our objective : to have a ``representative'' low-z early-type galaxy sample? --> sounds like a good idea.

I did a quick check of the FP relations using Meert catalog SerExp model fits. After matching
``two-component galaxies'' with H11 catalog, I selected galaxies with
probaE>0.7 & probaEll/probaE > 0.5 (about 2000 galaxies with z<0.05).
I narrowed the redshift (0.04 < z < 0.05) range so as not to worry about cosmology affecting distribution (~900 galaxies).

- L vs sigma (Faber-Jackson) : bulge and disk components have essentially
  identical relations ---> seperately virialised??
- surface brightness vs size (Kormendy) : 
- size vs sigma : 

Jenny thinks other than the Faber-Jackson, these scaling relations don't tell
us much other than that disk are the more extended component as supposed to be.
The only thing that could have helped us in arguing that we _need_ a second
component was **flattening**. We did a check of the distribution of
ellipcities of bulge and disk, and they are very similar instead of disk being
much flatter.

Of ~150 galaxies in the subsample that are in Meert catalog, 30 per cent have
``two-component'' flag set.

Plans
    - On the fitter side,
        * I need to get the two component fit working on the ~260 subsample,
          and compare with Meert results.
    - On the sample side,
        * I need to define the main and supplement (nearby and bright) samples
          and see their distribution in key parameters
    - Instead of arguing that we absolutely need a second component, we can
      investigate those which _can_ be fit by two component model, and those
      not.
        * Would stacking of images reveal differences in brightness profile?


2014-09-02 Tue meeting with Jenny

1. size of cutout images: 5-10 Re of single Sersic model?
    - not ideal since I would have to do an initial run of sersic on every
      galaxy just to determine the cutout sizes of real science images
2. masking -- since we cover large areas,
    - masking may not impact the model parameters at all..
    - we could find the location, and do aggressive masking



- currnet cutout size is not large enough for many galaxies
- check those that does not reach background (large) in Meert catalog
- freesky option

andrea kra

BGCs SDSS out to 100kpc halo/stellar doesn't flatten (26 mag/arcsec^2)

2014-07-07

## matching Jenny's ellipse results



2014-06-18
skype meeting w/ JG and CL

- cutoff keyword
- two examples of good/bad exponential component
  + size of inner/outer (use DVC)
  + half of the flux in outer part
  + ivar images: downweighted in the center flux^2, flux^1.5
  + higher sersic index in the outskirts

kinmetetry
- gave three example galaxies J1006, J2356, J0937 to Jenny with ivar images

2014-06-04
meeting w/ JG
- the exponential is not fitting the outskirts well...

2014-06-03
```
I've added (to the master copy) a python script called make_condor_shell.py. This is an update of an older perl script, so it may be buggy (I don't have condor here, so can't check it). You'll probably want to change the environment variables (you'll need to pass in the profile names). In any case, you need 3 parts to make the condor submit work.

1) shell file that actually runs the IDL code. This will run a specified number of galaxies. Don't make this number too small, or the overhead from condor won't make this worth it (usually you need jobs longer than 5 minutes). This needs the proper environment setup (starting with the shell). I think I've done that, but there's a chance the ${HOME} won't work, so put your own directory here and edit anything else you think you need to source (maybe the photoop tools from sdss?) At least, you'll need the right things in IDL_PATH. This can be finicky, and I remember the error messages being less than helpful.

2) condor submit file for each shell script. This is what condor needs to run the job. Note that you'll need folders ./condor/ and /tmp/condor/ for the error/output/log files. The log files have to live on your machine's local scratch space, not the NFS file system, for reasons I don't understand, but the jobs don't work otherwise. Along the same lines, you need to have the data somewhere all the machines can get to it. There are larger scratch spaces in the department (I think I used /peyton/scr/chimera), so ask. You can also read/write to the NFS. Check that this works by submitting a single condor job ($ condor_submit submit.0001 ). This is a good test that the environment setup is good, as well. If you use the /cutoff option, the galaxies should run a lot faster for debugging.

3) dag_submit. This allows you to submit lots of jobs at once. I usually ran ~50 at a time, but beware you may run out of IDL licenses. The maxjobs sets the total number of concurrent jobs.

To submit, just run
$condor_submit_dag dag_submit

This is probably not the only way to do this, but condor is fussy, so I ended up sticking with something that worked.
```

This makes my filetree quite messy. One easy alternative is using jug, and
indeed this works well, but the problem with using jug on condor is I cannot
seem to ask for a pool of cpus for vanilla universe -- understandably since
cpus are distributed among different machines. Gave up this option, but jug is
still useful and easy on my machine


2014-05-15

Meeting with Jenny
- We are not trying to come up with the "best model" of the surface brightness.
Instead, we want to test physically motivated core+outskirt decomposition using the
FP. Multi-component fits are always subject to degeneracies, and what we try here
is to use the independent information of velocity dispersion and the FP to minimize
(break?) this degeneracies.
- ellipticity & PA profile: It will take quite a while if I write my own code. Two
options available
    1. pyraf
    2. kinemetry (IDL) by Krajnovic http://davor.krajnovic.org/idl/#kinemetry
- we still want to test how robust our results are if we release the index from 4.



2014-05-14
started a separate note on model selection

2014-05-11

comments on the summary to J&C

Jenny

1. Some random questions -- in the fit images, are masked things still in
   there?  I ask because there are some that are mostly truncated.
2. Use ellipticity and PA as measures of goodness-of-fit
3. Use the FP to reduce degeneracies

Claire

1. How big are the chi-square values? Are they just driven by the shear number
   of pixles in these images?
2. Do an F-test between single deV and deV+Exp model. Technically, Keep in
   mind, the F-test is designed for _linear_ fits (not Sersic profiles), with
   _uncorrelated, Gaussian errors_, so all of the model testing we do is kind
   of cheating.

Steve sent me php codes that could be used to run image classifying program on
the server.

F-test : sensitive to non-normality?
Levene's test Bartlett's test
Brown-Forsythe test

2014-05-08

Liddle 2007 Information criteria for astrophysical model selection

Two schools in model selection

1. Beyesian inference
    - Bayesian Information Criterion (BIC, Schwarz 1978)
        + assumes the data points are independent and identically distributed
        + as an easier way to estimate Bayesian evidence?
    - DIC
2. Information-theoretic methods 
    - Akike Information Criterion (AIC) Takeuchi Information Criterion (TIC)

2014-05-07

model selection

* what is usually done (Huang+13)
    - examine residual image
        + finds obvious bad models insensitive to subtle differences possibly
        + subjective?
    - some goodness-of-fit quantity
        + ease of manual labor for a large sample objective uses integrated
        + information
    - compare 1D surface brightness profile
        + relatively sensitive to mild differences
* Huang+13
    - best model: minimum number of components with reasonable, robust
    - parameters that describes visibly distinct structure Excess Variation
    - Index (Hoyos+11)



2014-02-07

## Experiment following Zhu+10
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
