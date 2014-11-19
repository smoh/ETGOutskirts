## Resample PA and ellipticity curves from IRAF ellipse

For 32 galaxies, I test how PA/ellipticity curves from IRAF ellipse change by resampling each galaxy image 100 times. Each galaxy has a subdirectory containing the resampled images [0-99].fits. I run runall.cl script which runs all runJXXXX.cl scripts (one for each galaxy). In each subdirectory, each resampled image has corresponding ellipse output files in stsdas table format [0-99].out and this converted to ascii by tprint [0-99].dat.

Note that this failed for J0854 because my masking routine masked out the center of the galaxy. I restart manually after commenting out lines in runall.cl.

- resample.py : resample images and arange them in each subdirectory and write a cl script that will run ellipse task on all images
- makeplots.py : make plots of the result

20141119
fits images and *.out binary table deleted to save space 
