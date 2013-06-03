FUNCTION cauchyprior,r,x0,gamma

	out=x0+gamma*tan(!DPI*(r-0.5d))

	return,out
END
