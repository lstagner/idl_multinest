FUNCTION example_log_likelihood,params

    x=[-.35d0,-.4d0,-.2d0,.45d0,.1d0]
    y=[.2d0,-.4d0,.15d0,.1d0,-.15d0]
    sig=[.01d0,.01d0,.03d0,.05d0,.02d0]
    A=[1.0d0,.5d0,.8d0,.6d0,.5d0]
                
    d=A*exp(-((x-params[0])^2.0d + (y-params[1])^2.0d)/(2.0d * sig^2.0d))

    w=where(d eq 0.0d or d eq !values.f_nan,nw)

    if nw ne 0 then d[w]=(MACHAR(/double)).xmin

    tmp=total(d)

	return,alog(tmp)
END
                                    