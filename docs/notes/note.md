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
