FUNCTION gaussianprior,r,mu,sigma

	if (r le 1.0d-16) or (1.0-r) le 1.0d-16 then out=-1.0d32 else begin
		out=mu+sigma*sqrt(2.0d)*inverfc(2.0d*(1.0d-r))
	endelse

	return,out
END
