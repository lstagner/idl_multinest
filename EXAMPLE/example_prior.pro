FUNCTION example_prior,params

	tmp=params
    tmp[0,*] = -1.0d + 2.0d * (tmp[0,*])
    a= -( 1.0d - (tmp[0,*])^2.0d )^.5d
    b=  ( 1.0d - (tmp[0,*])^2.0d )^.5d
    tmp[1,*]= a[0] + (b[0]-a[0])*tmp[1,*]
    return,tmp
                                                	
END