PRO test_clustering,NUM=NUM

    DefSysV, '!RNG', Obj_New('RandomNumberGenerator')
	if not keyword_set(NUM) then NUM=10
	
	for i=0,NUM-1 do begin
		if i eq 0 then begin
			r=!RNG->GetRandomNumbers(2,5,/NORMAL)
			r[0,*]=.2*r[0,*]+0.0*r[1,*]-.1
			r[1,*]=0.0*r[0,*]+.2*r[1,*]-.1 
			;print,transpose(r)
			;print,''
		endif else begin
        	d=!RNG->GetRandomNumbers(2,30,/NORMAL)
        	d[0,*]=.2*d[0,*]+.5
        	d[1,*]=.2*d[1,*]+.5
        	;print,transpose(d)
        	r=[[r],[d]]
        endelse
    endfor 
    
    clusters=cluster_data(r)
    
   	num_clusters=n_elements(clusters)
   	
	window,0,xsize=500,ysize=500   	
   	loadct,0,/silent
   	plot, r[0,*],r[1,*],psym=1,xrange=[-2,2],yrange=[-2,2]
  	loadct,13,/silent
   	for i=0,num_clusters-1 do begin
   		mean=clusters[i].mean    
   		cov=clusters[i].cov
   		invcov=clusters[i].invcov
   		evals=clusters[i].evals
   		evecs=clusters[i].evecs
   		
   		d_ptr=clusters[i].data_ptr
   		d=*d_ptr    
   		k=dblarr(clusters[i].data_num)
  
   		for j=0L, clusters[i].data_num-1 do begin
   			tmp=d[*,j]-mean
   			k[j]=tmp##invcov##transpose(tmp)
   		endfor
 
   		kmax=max(k,nk)
   		print,d[*,nk]
   		for j=0,n_elements(evals)-1 do begin
   			evecs[*,j]=sqrt(evals[j])*evecs[*,j]
   		endfor
   		
		rd=dblarr(2,100000)
		for j=0L, n_elements(rd[0,*])-1 do begin
			rd[*,j] = sample_ellipsoid(cov,mean,kmax,expand=1.08,scaled_evecs=evecs)
	   	endfor
	   
	   	oplot,rd[0,*],rd[1,*],psym=1,color=50+i*50
	   	loadct,0,/silent
	   	oplot,r[0,*],r[1,*],psym=1
	   	print,ellipsoid_vol(clusters[i].detcov,1.08*kmax,n_elements(mean))
	   	
   	endfor
END
