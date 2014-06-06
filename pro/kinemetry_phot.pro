;------------------------------------------------------------------------------
;NAME:
;   KINEMETRY_PHOT
;PURPOSE:
;   Simplify the original KINEMETRY function for fitting surface brightness
;	images
;CALLING SEQUANCE:
;	KINEMETRY_PHOT, IMAGE, X0, Y0, ANGIN, QIN, OUTNAME, NTRM=10, /VERBOSE, /PLOT
;OUTPUT:
;	Saves a table of RAD, PA, Q, CF, ER_PA, ER_Q, ER_CF to a FITS table
;PARAMETERS:
;	IMAGE : file name of the image
;	X0, Y0 : center position, fixed during ellipse fitting
;	ANGIN, QIN : initial values for ANGIN, QIN
;	OUTNAME : output FITS name
;	NTRM : number of the Fourier terms
;	/VERBOSE : print verbose output
;	/PLOT : set plot keyword for kinemetry
;------------------------------------------------------------------------------

PRO KINEMETRY_PHOT, IMAGE, X0, Y0, ANGIN, QIN, OUTNAME, $
        VERBOSE=VERBOSE, NTRM=NTRM, PLOT=PLOT

ext = 0
img = mrdfits(image, ext, h)
s = size(img)

; 
; determine size and make dummy 1-D arrays
;
n=s[1]*s[2]
yy=REPLICATE(1,s[1])#(indgen(s[2]))
xx=(indgen(s[1]))#REPLICATE(1,s[2])
x=REFORM(xx, n)
y=REFORM(yy, n)
flux=REFORM(img, n)

;
; running kinemetry
;
KINEMETRY, x, y, flux, rad, pa, q, cf, NTRM=NTRM, /EVEN, verbose=VERBOSE, $
           X0=X0, Y0=Y0, XC=xc, YC=yc, IMG=img, $
           ER_CF=er_cf, ER_PA=er_pa, ER_q=er_q, ER_XC=er_xc, ER_YC=er_yc, $
           XELLIP=xellip, YELLIP=yellip, $
           /nogrid, PAQ=[ANGIN, QIN], /fixcen, PLOT=PLOT

; Save output
struct = {rad:rad, pa:pa, q:q, cf:cf, er_pa:er_pa, er_q:er_q, er_cf:er_cf}
mwrfits, struct, outname, /create

END
