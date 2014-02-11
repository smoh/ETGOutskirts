;+
; NAME: 
;          SGN
;
; PURPOSE: 
;          return the sign (+1 or -1) of a number
;
; CALLING SEQUENCE:
;
;          sign = sgn(num)
;
; INPUTS:
;
;          num: number
;
; OPTIONAL INPUTS:
;
;          NONE
;
; KEYWORD PARAMETERS:
; 
;          NONE
;
; OUTPUTS:
;
;         -1 or +1 as an INTEGER
;
; EXAMPLE:
;         IDL> print, sgn(-4.)
;          -1
;
; MODIFICATION HISTORY:
;         Apr 13 2010
;         Claire Lackner
;-

;returns the sign (-1 or +1) of a number, 0 is plus 1
FUNCTION sgn, num
neg = ( num lt 0 ) * (-1)
pos = (num gt 0 or num eq 0)
return, neg+pos
end
