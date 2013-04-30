FUNCTION ellipsoid_vol,detcov,k_fac,ndim,expand=expand

	if not keyword_set(expand) then expand=1.0D
	PI=3.14159265D0
	
	ellvol=sqrt( detcov*( (k_fac*expand)^(ndim) ) )
	
	if (ndim mod 2) eq 0 then begin
		i=2.0D
		while (i le ndim) do begin
			ellvol=ellvol*2.0D * PI / i
			i=i+2.0D
		endwhile
	endif else begin
		ellvol=2.0D * ellvol
		i=3.0D 
		while (i le ndim) do begin
			ellvol = ellvol * 2.0D * PI / i
			i=i+2.0D
		endwhile
	endelse
	
	return, ellvol
END