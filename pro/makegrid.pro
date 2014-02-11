;+
; NAME: 
;          MAKEGIRD
;
; PURPOSE: 
;          make two arrays of (float) grid indices given the desired dimensions
;
; CATEGORY:
;          
;          array
;
; CALLING SEQUENCE:
;
;          makegrid, xlen, ylen, xgrid, ygrid
;
; INPUTS:
;
;          xlen, ylen:    integers giving size of xdimension of grid
;          xgrid, ygrid:  empty variables in which to store grided arrays
;
;
; KEYWORD PARAMETERS:
;     x_size: spacing of points in x, so the total number of points is xlen
;     y_size: spacing of points in y, so the total number of points is ylen
;     x_start, y_start: coordinates of the first point, default is (0,0)
;
; EXAMPLE:
;
;         to draw a circle on a blank array
;         IDL> xl = 5
;         IDL> yl = 5
;         IDL> makegrid( xl, xl, xgrid, ygrid )
;         IDL> print, ygrid
;          0.00000      0.00000      0.00000      0.00000      0.00000
;          1.00000      1.00000      1.00000      1.00000      1.00000
;          2.00000      2.00000      2.00000      2.00000      2.00000
;          3.00000      3.00000      3.00000      3.00000      3.00000
;          4.00000      4.00000      4.00000      4.00000      4.00000
;         IDL> print, xgrid
;          0.0000000       1.0000000       2.0000000       3.0000000       4.0000000
;          0.0000000       1.0000000       2.0000000       3.0000000       4.0000000
;          0.0000000       1.0000000       2.0000000       3.0000000       4.0000000
;          0.0000000       1.0000000       2.0000000       3.0000000       4.0000000
;          0.0000000       1.0000000       2.0000000       3.0000000       4.0000000
;
; MODIFICATION HISTORY:
;         Claire Lackner, Feb 4 2010
;-
PRO makegrid, xlen, ylen, xgrid, ygrid, $
              x_size=x_size, y_size=y_size, $
              x_start=x_start, y_start=y_start

row= dindgen(xlen)
col = dindgen(ylen)

if( keyword_set(x_size) ) then row*=x_size
if( keyword_set(y_size) ) then col*=y_size

if(keyword_set(x_start) ) then row += x_start
if(keyword_set(y_start) ) then col += y_start

xgrid = row#(col*0+1)
ygrid = (row*0+1)#col

END
