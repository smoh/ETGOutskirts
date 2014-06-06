
massive.dat -- all galaxies in massive survey

massive_SDSS.fits: matching with NSA
- 70/106 galaxies matched within 5 arcsec separation
- maximum separation ~1.8 arcsec

download SDSS and NSA images
$downloadNSA massive_SDSS.fits data/nsa NSA
$downloadNSA massive_SDSS.fits data/sdss_psf_meta psField
$downloadNSA massive_SDSS.fits data/sdss_field r

prepare SDSS images for bdfitter
$python ../subsample_test/prepare_sdss.py massive_SDSS.fits data
took about 2hrs

