;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
; Jul 19 2010
;
; secant.pro
;
; solves an equation using the secant method
;;;;;;;;;;;;;;;;;;;;;;;;
;+
;NAME: SECANT
;
;CALLING: secant(func, xn_1, xn, precision=precision, maxiter=maxiter)
;
;INPUT: func = function to be zerod
;       xn_1 = initial guess first value
;       xn   = initial guess second value
;
;KEYWORDS: precision = stop when |x_j - x_(j-1)| < precision
;                      (default = 1.0e-5)
;          maxiter = maximum number of steps (default = 100)
;
;OUTPUT: x where func(x)=0 within precision
;-

FUNCTION testfunc, x
return, cos(x)-x^3
END

FUNCTION secant, func, xn_1, xn, precision=precision, maxiter=maxiter

;set keyword default values
if not keyword_set(maxiter) then maxiter=100
if not keyword_set(precision) then precision=1.0e-8

for n=1,maxiter-1 do begin
;    print, xn
    d=(xn-xn_1)/(call_function(func,xn) - $
                 call_function(func,xn_1)) * $
      call_function(func,xn)

    bad=where(abs(d) ge precision, nbad)

    if nbad eq 0 then return, xn
    xn_1 = xn
    xn = xn - d
endfor

print, "secant method: maxiter reached"
return,xn


END
