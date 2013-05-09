FUNCTION example_prior,params
	n_samps=n_elements(params)
	output=ptrarr(n_samps)
	for i=0,n_samps-1 do begin
		tmp=*params[i]
	    tmp[0] = -1.0d + 2.0d * tmp[0]
	    a= -( 1.0d - (tmp[0])^2.0d )^.5d
	    b=  ( 1.0d - (tmp[0])^2.0d )^.5d
	    tmp[1]= a[0] + (b[0]-a[0])*tmp[1]
		output[i]=ptr_new(tmp)
	endfor
    return,output
                                                	
END
