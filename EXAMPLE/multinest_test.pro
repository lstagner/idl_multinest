;;WHEN RUN THIS EXAMPLE (MULTINEST_TEST) SHOWS THE CLUSTERING PROCESS FOR THE EXAMPLE FUNCTION

;;NEEDED FUNCTIONS
FUNCTION log_plus , a , b

	if a gt b then begin
		tmp = double(a) + alog(1.0+exp(double(b)-double(a)))
	endif else begin
		tmp = double(b) + alog(1.0+exp(double(a)-double(b)))
	endelse
	
	return, tmp
END

FUNCTION in_prior,point
	w=where(point gt 1.0D or point lt 0.0D,nw)
	if nw ne 0 then return, 0 else return,1
END

FUNCTION prior, params
 
 		tmp=params
		tmp[0,*] = -1.0d + 2.0d * (tmp[0,*])
		a= -( 1.0d - (tmp[0,*])^2.0d )^.5d
		b=  ( 1.0d - (tmp[0,*])^2.0d )^.5d	
		tmp[1,*]= a[0] + (b[0]-a[0])*tmp[1,*]
		return,tmp
END

FUNCTION log_likelihood,params
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

FUNCTION remove_value,index,array
                       
    num=n_elements(array)
    
    if index eq 0 then new_array=array[1:num-1]
    if index eq num-1 then new_array=array[0:num-2]
                           
    if index ne 0 and index ne num-1 then begin
    	new_array=[array[0:index-1],array[index+1:num-1]]
    endif
    
    return,new_array
	
END

PRO multinest_test,NUM_PARAMS, NUM=NUM,SAMPLE_NUM=SAMPLE_NUM,EXPAND=EXPAND,PLOT=PLOT
	
	DefSysV, '!RNG', Obj_New('RandomNumberGenerator')
	if not keyword_set(NUM) then NUM=1000
	if not keyword_set(SAMPLE_NUM) then SAMPLE_NUM=1
	if keyword_set(PLOT) then window,0,xsize=600,ysize=600

	;;GET INITIAL POP
	live_points = !RNG -> GetRandomNumbers(NUM_PARAMS,NUM,/double)
	
	samples=prior(live_points)
	
	;;SET INITIAL LOG LIKELIHOODS
	likelihoods=dblarr(NUM)
	for i=0L, NUM-1 do begin
		likelihoods[i]=log_likelihood(samples[*,i])
	endfor
	
	likelihood_min=MIN(likelihoods,likelihood_min_loc)
	H = 0.0

	logZ=double(-10.0^(30.0))
	k=0.0
	logw=alog(1.0-exp(-1.0/double(NUM)))
	
	while (k le 2.0*H*double(NUM)) do begin		
		;;CLUSTER DATA
		result=cluster_data(live_points)
		clusters=n_elements(result)
		;;FIGURE OUT HOW MANY POINTS IN EACH CLUSTER
		ellnum=dblarr(clusters)
		for i=0,clusters-1 do begin
			ellnum[i]=result[i].data_num
		endfor
		
		;;GET SAMPLE_NUM NUMBER OF SAMPLES FROM THE ELLIPSOIDS		
		totnum=total(ellnum)
		cnt=0L
		while cnt le SAMPLE_NUM-1 do begin
			;;PICK A RANDOM CLUSTER AND SAMPLE IT
			ellprob=total(double(ellnum)/double(totnum),/cumulative)
			r=!RNG->GetRandomDigits(1)
			w=in_ellipsoids(live_points[*,(r mod NUM)],result, expand=EXPAND, /location)
			
			if w[0] eq -1 then continue
			if n_elements(w) gt 1 then w[0]=w[r mod n_elements(w)]
			cov=result[w[0]].cov
			mean=result[w[0]].mean
			sevecs=result[w[0]].sevecs                        
			k_fac_max=result[w[0]].kmax
			;;SAMPLE CLUSTER UNTIL VALID POINT IS FOUND
			while 1 do begin
				samp=sample_ellipsoid(cov,mean,k_fac_max,expand=EXPAND,scaled_evecs=sevecs)
				prior_check=in_prior(samp)
				if prior_check eq 0 then continue
				tsamp=prior(transpose(samp))
				samp_likelihood=log_likelihood(transpose(tsamp))
				if (samp_likelihood gt likelihood_min) then begin
	        		break
				endif
			endwhile
			
			;;CHECK TO SEE IF SAMPLE IS IN ANY OTHER CLUSTER AND RANDOMLY ACCEPT IT BASED
			;;ON HOW MANY CLUSTERS IT IS IN
			samp_loc=in_ellipsoids(transpose(samp),result,expand=EXPAND)
			r=!RNG->GetRandomNumbers(1,/double)
			if r le (1.0d/samp_loc) then begin
				;;UPDATE SAMPLING DATA 	    	                                                                                    	    				
				samples=[[samples],[tsamp]]
				live_points=[[live_points],[transpose(samp)]]
				likelihoods=[likelihoods,samp_likelihood]
				logwt=double(logw+likelihood_min)
				logZnew=log_plus(logZ,logwt) 
				H = exp(logwt-logZnew)*likelihood_min + exp(logZ-logZnew)*(H+logZ)-logZnew
				logZ = logZnew
				logw = logw-1.0/(double(NUM))
				
				w=where(likelihoods ne likelihood_min,w_num,complement=nw,ncomplement=nw_num)
				if nw_num gt 1 then w=[w,nw[1:*]]
			    likelihoods=likelihoods[w]
				live_points=live_points[*,w]
				                                
				likelihood_min=MIN(likelihoods,likelihood_min_loc)
		
				cnt=cnt+1

			endif
		endwhile
		

		k=k+cnt
		if k mod 100 eq 0 then print,100.0d*k/(2.0d * H * double(NUM))," % Completed "

		;;PLOT OUT CURRENT LIVE POINTS AND SHOW HOW THEY ARE CLUSTERED
		if keyword_set(PLOT) then begin
			new_samp=dblarr(2,1000)
			for m=0,clusters-1 do begin
		   		evecs=result[m].sevecs                 
				;;DRAWS ELLIPSOID OUTLINE                
		 		for j=0L,999 do begin 
		 			tmp=sample_ellipsoid(result[m].cov,result[m].mean,result[m].kmax,expand=EXPAND,scaled_evecs=evecs,/surface)
		 			new_samp(*,j)=tmp
		 		endfor
		 		loadct,0,/silent
			 	if m eq 0 then plot,new_samp(0,*),new_samp(1,*),psym=3,xrange=[0,1],yrange=[0,1] else $
				 	oplot,new_samp(0,*),new_samp(1,*),psym=3
				;;DRAWS OUT POINTS IN ELLIPSOID
			 	loadct,13,/silent
			 	d_ptr=result[m].data_ptr
			 	d=*d_ptr
				oplot,d[0,*],d[1,*],psym=3,color=50*m+50 
			endfor
		endif
		for i=0,clusters-1 do ptr_free,result[i].data_ptr

	endwhile
	
	;;UPDATE SAMPLING DATA FOR REMAINING POINTS
	logw=-double(k)/double(NUM)-alog(double(NUM))
	for i=0, NUM-1 do begin
		logwt=double(logw+likelihoods[i])		
		logZnew=log_plus(logZ,logwt)
		H = exp(logwt-logZnew)*likelihoods[i] + exp(logZ-logZnew)*(H+logZ)-logZnew
		logZ=logZnew	
	endfor
	
	print,"H = " , H 
	print,"Log(Z) = ", logZ , " +- ", sqrt(H/(1.0*NUM))
	if keyword_set(PLOT) then begin
		loadct,0,/silent
		plot,samples[0,*],samples[1,*],psym=3,xrange=[-1,1],yrange=[-1,1]
	endif
	return, samples
	GET_OUT:
END

  
