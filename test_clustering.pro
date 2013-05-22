PRO test_clustering,NUM=NUM,time=time,k_medoids=k_medoids

    DefSysV, '!RNG', Obj_New('RandomNumberGenerator')
	if not keyword_set(NUM) then NUM=10
	if not keyword_set(time) then begin
		for i=0,NUM-1 do begin
			if i eq 0 then begin
				r=!RNG->GetRandomNumbers(2,30,/NORMAL)
				rr=!RNG->GetRandomNumbers(2,/double)
				rr=-10.0+20.*rr
				r[0,*]=r[0,*]+0.0*r[1,*]+rr[0]
				r[1,*]=0.0*r[0,*]+r[1,*]-rr[1]
				;print,transpose(r)
				;print,''
			endif else begin
	   	     	d=!RNG->GetRandomNumbers(2,30,/NORMAL)
				rr=!RNG->GetRandomNumbers(2,/double)
				rr=-10.0+20.*rr
   		     	d[0,*]=d[0,*]+rr[0]
        		d[1,*]=d[1,*]+rr[1]
        		;print,transpose(d)
        		r=[[r],[d]]
       		 endelse
   		endfor
	    clusters=cluster_data(r,k_medoids=k_medoids)
	   	plot_2d_clusters,clusters,x_range=[-10,10],y_range=[-10,10]
	endif else begin
		nclusters=[1,2,5]
		npoints=[10,100,1000,10000]
		times=dblarr(n_elements(npoints),n_elements(nclusters))
		for m=0,9 do begin
		for j=0,n_elements(npoints)-1 do begin
			for k=0,n_elements(nclusters)-1 do begin
				for i=0,nclusters[k]-1 do begin
					if i eq 0 then begin
						r=!RNG->GetRandomNumbers(2,long(npoints[j]/nclusters[k]),/NORMAL)
						rr=!RNG->GetRandomNumbers(2,/double)
						rr=-10.0+20.*rr
						r[0,*]=r[0,*]+0.0*r[1,*]+rr[0]
						r[1,*]=0.0*r[0,*]+r[1,*]-rr[1]
						;print,transpose(r)
						;print,''
					endif else begin
		   	    	 	d=!RNG->GetRandomNumbers(2,long(npoints[j]/nclusters[k]),/NORMAL)
						rr=!RNG->GetRandomNumbers(2,/double)
						rr=-10.0+20.*rr
   			     		d[0,*]=d[0,*]+rr[0]
   		     			d[1,*]=d[1,*]+rr[1]
  	    	  			;print,transpose(d)
        				r=[[r],[d]]
       				endelse
   				endfor
				t=systime(1)
	    		clusters=cluster_data(r,k_medoids=k_medoids)
				times[j,k]=times[j,k]+systime(1)-t
			endfor
		endfor
		endfor		 
    	print,times/10.0
	endelse
END
