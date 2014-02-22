profiles={SER:8}   ;this order is important! (=order in model image FITS)
; Be sure to include trailing '/' in directories
FIT_SAMPLE, 'input_SER.fits', 0, 51, './out_SER/', './data/',  $
    /residuals, profiles=profiles

