FUNCTION logprior,r,x1,x2

	if r le 0.0d then out = -1.0d32 else begin
		lx1=alog(x1) 
		lx2=alog(x2)
		out=exp(lx1+r*(lx2-lx1))
	endelse
	
	return, out
END
