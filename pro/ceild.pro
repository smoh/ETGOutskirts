;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Claire Lackner
; Sept 15, 2009
; ceilings to nearest double with 10^digits as last significant digit 
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;rounds up to the nearest double with 10^digits as last significant
;digit 
function ceild, num, digits
pow=10.^digits
n=ceil(num/pow)*pow
return, n
end

