# Subsample test

2014-03-10
excluded J101933 (ring) and J160339 (close companion)

Fitting
-------

single_DVC/
single De Vaucouleurs profile

single_SER/
single Sersic profile
```
$../fit_sample --res --debug input_SER.fits 0 48 out ../data > log.txt
$ python ../plot.py ../SampleZMprobaEllSub_visual.fits out/RAWFIT00000.00048.fits DVC ../data out/models plots
$ python ../make_webimages.py out/RAWFIT00000.00048.fits ./plots summary.html
```


