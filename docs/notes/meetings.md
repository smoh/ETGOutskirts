
2014-02-10 Skype w/ Claire on setting things up

- ivar: inside parent images
- take a closer look at `make_input.py (I could write my own version of this)`, 
  `fit_sample.pro`, `sersic.pro`
- she will be working on giving upper/lower limit to parameters (rather than just fixed). This is not high priority at the moment.
- should take a look at condor system at Peyton
- the main dependency of the whole code is idlutils.. but maybe others, too.
- Let's start with a small sample of 10-50 that covers a range of redshift, sigma, $M_r$, $R_e$

2014-03-06 Jenny

In any case, we need to infer size from dispersion. Try two things:
1. sigma --(Faber-Jackson relation)--> L --(Kormendy relation)--> size
2. define our own size-velocity dispersion relation using SDSS deVac model mag (see van der Wel et al. 2008 -- large scatter): median / some smaller size constraint for compact objects


