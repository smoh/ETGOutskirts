header

data/
  sdss_field/
  sdss_psf_meta/
  nsa/
    pimages/

- Make stamp images from field images
$etg_make_stamps sample.fits r --outdir data > cutout.log
output:
data/
  sdss_ivar/
  cutout/
  images/
  ivar/
  mask/
  psf/
  fig/stamps

- Do a single Sersic fitting for additional masking
in testser/
$python gen_input.py --> input_ser.fits
my.exec, my.job copied from previous runs (10 gals / queue)
condor job took at most ~1h 30min

```
J083121.03+461314.5:            9 of           10
image: data/images/J083121.03+461314.5.fits
 ivar: data/ivar/J083121.03+461314.5.fits
MRDFITS: Image array (769,817)  Type=Real*4
MRDFITS: Image array (769,817)  Type=Real*4
MRDFITS: Image array (51,51)  Type=Real*4
   cropped size:         769x         817
doing profile SER
using given IC
       1.0000000       10.000000       4.0000000      0.80000000       0.0000000
       358.97600       407.80800       0.0000000       0.0000000
Iter      1   CHI-SQUARE =       5457586.5          DOF = 628266
    P(0) =              1.00000
    P(1) =              10.0000
    P(2) =              4.00000
    P(3) =             0.800000
    P(4) =              0.00000
    P(5) =              358.976
    P(6) =              407.808
    P(7) =              0.00000
    P(8) =              0.00000
Iter      2   CHI-SQUARE =       5371457.8          DOF = 628266
    P(0) =          1.11022E-16
    P(1) =              13.2797
    P(2) =              4.83759
    P(3) =             0.827692
    P(4) =              0.00000
    P(5) =              359.022
    P(6) =              407.828
    P(7) =           0.00544119
    P(8) =              0.00000
Iter      2   CHI-SQUARE =       5371457.8          DOF = 628266
    P(0) =          1.11022E-16
    P(1) =              13.2797
    P(2) =              4.83759
    P(3) =             0.827692
    P(4) =              0.00000
    P(5) =              359.022
    P(6) =              407.828
    P(7) =           0.00544119
    P(8) =              0.00000
```

- Do additional masking
$etg_admask testser testser/RAWFIT00000.00268.fits data testser/models \
data/admask data/fig/admask > admask.log
--> data/admask, data/fig/admask

(This etg_admask is copied and modified from test2/flag.py)

- ExamineStamps.ipynb 
--> Updated sample.fits with fracarea, cutout_cx, cutout_cy, fracRe, frac2Re

## Fitting models

All fittings are in fit/

- Ser/ -- single Sersic
Use best-fit params of test Sersic fit as initial values
Combined RAWFIT files --> fit_Ser.fits

- deV/ -- single deV
Use best-fit params of test Sersic fit as initial values, fixing n=4
Combined RAWFIT files --> fit_deV.fits

- deVExp -- deV + exp model

- n2Exp -- n=2 + exp model

## yang group catalog

- `imodel[A,B,C]_combined` : `imodelA_1` matched with `SDSS7` by exact value galaxyid. groupN column is
    generated from internal matching of groupid

- what are group id = 0?
  - modelA has no groupid=0 entry
- how different is model A, B, C?
- Is yang modelC the same as the catalog Jenny gave me?
  - No.
- `sample_imodelC.csv` : match db.csv ra, dec with `imodelC_combined` best
  match symmetric, join all from db.csv.

## stellar age catalogs


