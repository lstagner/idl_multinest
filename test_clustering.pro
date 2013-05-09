PRO test_clustering,NUM=NUM

    DefSysV, '!RNG', Obj_New('RandomNumberGenerator')
	if not keyword_set(NUM) then NUM=10
	
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
    
    clusters=cluster_data(r)
   	plot_2d_clusters,clusters,x_range=[-10,10],y_range=[-10,10]
END
