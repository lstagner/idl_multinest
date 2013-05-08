FUNCTION example_prior,params

	tmp=params
	n_samps=n_elements(params)
	for i=0,n_samps-1 do begin
	    (*tmp[i])[0] = -1.0d + 2.0d * *(tmp[i])[0])
	    a= -( 1.0d - ((*tmp[i])[0])^2.0d )^.5d
	    b=  ( 1.0d - ((*tmp[i])[0])^2.0d )^.5d
	    (*tmp[i])[1]= a[0] + (b[0]-a[0])*(*tmp[i])[1]]
	endfor
    return,tmp
                                                	
END
