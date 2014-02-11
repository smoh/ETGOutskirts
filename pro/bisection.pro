;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
; Jul 19 2010
;
; secant.pro
;
; solves an equation using the bisection method
;;;;;;;;;;;;;;;;;;;;;;;;
;+
;NAME: BISECTION
;
;CALLING: bisection(func, xmin, xmax, precision=precision, maxiter=maxiter)
;
;INPUT: func = function to be zerod
;       xmin = initial guess first value
;       xmax   = initial guess second value
;
;KEYWORDS: precision = stop when |xmax-xmin| < precision
;                      (default = 1.0e-5)
;          maxiter = maximum number of steps (default = 100)
;
;OUTPUT: x where func(x)=0 within precision
;-

FUNCTION testfunc, x
return, cos(x)-x^3
END

FUNCTION bisection, func, xmin, xmax, precision=precision, maxiter=maxiter

;set keyword default values
if not keyword_set(maxiter) then maxiter=100
if not keyword_set(precision) then precision=1.0e-5


for n=1,maxiter do begin
    xmid = (xmax+xmin)*0.5
    bad = where(xmax-xmin ge precision, nbad)
    if( nbad eq 0 ) then return, xmid
    func_mid = call_function(func,xmid)
;    print, xmid, func_mid
    move_down = where(func_mid*call_function(func,xmin) lt 0.0, nmd)
    move_up = where(func_mid*call_function(func,xmax) lt 0.0, nmu)
    if(nmd gt 0 ) then xmax[move_down] = xmid[move_down]
    if(nmu gt 0) then xmin[move_up] = xmid[move_up]
    if(nmu eq 0 and nmd eq 0) then return, xmid
endfor

print, "bisection method: maxiter reached"
return, xmid


END
