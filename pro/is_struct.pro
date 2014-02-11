;is_struct.pro

;Returns true is passed parameter is a struct or array of structs
;2008-02-11
;Fergal Mullally

function is_struct, obj
	if size(obj, /type) eq 8 then begin
		return, 1
	endif else begin
		return, 0
	endelse
    end
