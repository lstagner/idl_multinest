FUNCTION lognormalprior,r,a,sigma

	b=sigma*sigma + sigma*sqrt(2.0d)*inverfc(2.0d * r)
	out=a*exp(b)

	return,out
END
