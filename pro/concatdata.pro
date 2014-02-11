;____________________________________
;
; Claire Lackner
;  Nov, 5, 2009
;
;  concatdata.pro
;
;  IDL
; reads the output tables from the fitter and puts them into a single
; file, arguments are fold of input, wildcard string of data files,
; name of output file
;   
;____________________________________


PRO concatdata, folder, infiles, outfile, nosort=nosort, $
                maxdepth=maxdepth, wholename=wholename

;usage statement
if( n_params() ne 3 ) then message, $
  'Usage:concatdata, folder, infiles, outputfile'

if (keyword_set(maxdepth)) then md=maxdepth else md=1
files = unixfind(folder, infiles, maxdepth=md)

if( not keyword_set(nosort) ) then begin
    filenames = files[sort(files)] ; automatically sort the files
;    print, 'got here'
endif else filenames = files
print, filenames

temp = mrdfits( filenames[0], 1)
data = temp

for i=1L, n_elements( filenames )-1 do begin
    temp = mrdfits( filenames[i], 1)
    data = [[data], temp]
endfor

mwrfits, data, outfile, /create

END


