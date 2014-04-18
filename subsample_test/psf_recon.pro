
psfield_file = 'data/sdss_psf_meta/psField-001035-2-0113.fit'
row = 484.33
col = 1643.35

pstruct = mrdfits(psfield_file, 3)
nrow_b=(pstruct.nrow_b)[0]
ncol_b=(pstruct.ncol_b)[0]
;assumes they are the same for each eigen
;so only use the 0 one
rnrow=(pstruct.rnrow)[0]
rncol=(pstruct.rncol)[0]

nb=nrow_b*ncol_b
coeffs=fltarr(nb)
ecoeff=fltarr(3)
cmat=pstruct.c

rcs=0.001
for i=0l, nb-1 do coeffs[i]=(row*rcs)^(i mod nrow_b) * (col*rcs)^(i/nrow_b)


for j=0, 2, 1 do begin
    print, j
    for i=0l, nb-1, 1 do begin
    print, j, i
        ecoeff[j]=ecoeff[j]+cmat(i/nrow_b,i mod nrow_b,j)*coeffs[i]
    endfor
endfor

psf = (pstruct.rrows)[*,0]*ecoeff[0]+$
      (pstruct.rrows)[*,1]*ecoeff[1]+$
      (pstruct.rrows)[*,2]*ecoeff[2]
end