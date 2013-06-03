FUNCTION sinprior,r,x1,x2,degree=degree

	if keyword_set(degree) then begin 
		x1*=!DPI/180.0d
		x2*=!DPI/180.0d
	endif

	cx1=cos(x1) & cx2=cos(x2)

	out=acos(cx1+r*(cx2-cx1))
	return,out
END
