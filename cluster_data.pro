;;cluster_data.pro: Performs X-means clustering
;;
;;Arguments:
;;	n x m matrix where n is the dimensionality of the data and m is the number of points
;;Returns:
;;	An array of structures. Each structure is a cluster and contains a pointer to the data and other useful 
;; 	information 

FUNCTION cluster_data,data_pts,ellipsoid=ellipsoid

	FORWARD_FUNCTION cluster_data

	data_info=size(data_pts)	
	if data_info[-2] eq 10 then begin
		ncol=n_elements(data_pts) ;;NUMBER OF DATA POINTS
		nrow=n_elements(*data_pts[0]) ;;NUMBER OF DIMENSIONS/PARAMS
		data=dblarr(nrow,ncol)
		for i=0,ncol-1 do data[*,i]=*data_pts[i]
	endif else begin
		ncol=data_info[2] ;;NUMBER OF DATA POINTS
		nrow=data_info[1] ;;NUMBER OF DIMENSIONS/PARAMS
		data=data_pts
	endelse

	if not keyword_set(ellipsoid) then begin
		mean1=dblarr(nrow)
		for i=0,nrow-1 do begin
			mean1[i]=mean(data[i,*])
		endfor
		;;DETERMINE COVARAINCE MATRIX
		cov1=correlate(data,/covariance,/double)
		;;DETERMINE INVERSE COVARIANCE MATRIX
		invcov1=la_invert(cov1,/double)

		;;DETERMINE EIGENVALUES AND EIGENVECTORS OF COVARIANCE MATRIX
		evals1=eigenql(cov1,/double,eigenvectors=evecs1) > (MACHAR(/double)).xmin
		;;DETERMINE THE DETERMINANT OF THE COVARIANCE MATRIX
		detcov1=product(evals1)
		;;SCALE THE EIGENVECTORS BY THE SQRT OF RESPECTIVE EIGENVALUES
		sevecs1=evecs1
		for i=0,n_elements(evals1)-1 do begin
			sevecs1[*,i]=sqrt(evals1[i])*evecs1[*,i]
		endfor
		
		;;DETERMINE THE K FACTOR
		a=0
		kfac1=dblarr(ncol)
		for i=0,ncol-1 do begin
			kfac1[i]=transpose(data[*,i]-mean1)#invcov1#(data[*,i]-mean1)
			a=a+kfac1[i]
		endfor
		kmax1=max(kfac1)
		;;DETERMINE ELLIPSOID VOLUME
		ellvol1=ellipsoid_vol(detcov1,max(kfac1),nrow)
	endif else begin
		cov1=ellipsoid.cov
		detcov1=ellipsoid.detcov
		evals1=ellipsoid.evals
		evecs1=ellipsoid.evecs
		sevecs1=ellipsoid.sevecs
		kmax1=ellipsoid.kmax
		ellvol1=ellipsoid.ellvol
		invcov1=ellipsoid.invcov
		a=ellipsoid.ksum
		mean1=ellipsoid.mean
	endelse
	
	;;DETERMINE BAYESIAN INFORMATION CRITERION. TAKEN FROM MULTINEST BY F. FEROZ
	BIC1=nrow*ncol*alog(2.0*3.14159265)+ncol*alog(detcov1)+a+.5*nrow*(nrow+3.0)*alog(ncol)	
	
	;;TRY TWO CLUSTERS
;	w2=CLUST_WTS(data,N_CLUSTERS=2,N_ITERATIONS=100)
	w2=k_medoids(data,n_clusters=2,n_iterations=100)
	r2=CLUSTER(data,w2,N_CLUSTERS=2)
	clusters=r2[UNIQ(r2,sort(r2))]
	nwh1=0
	nwh2=0
	if n_elements(clusters) gt 1 then begin
		wh1=where(clusters[0] eq r2,nwh1)
		wh2=where(clusters[1] eq r2,nwh2)
	endif 
	
	if nwh1 lt nrow+1 or nwh2 lt nrow+1 then begin
		BIC2=0
		ellvol21=ellvol1
		ellvol22=ellvol1 
	endif else begin
		;;CONSTRUCT ELLIPSOIDS FOR EACH CLUSTER
		d1=data[*,wh1]
		d2=data[*,wh2]
		d1_info=size(d1)
		d2_info=size(d2)
		ncol1=d1_info[2]
		nrow1=d1_info[1]
		ncol2=d2_info[2]
		nrow2=d2_info[1]
		mean21=dblarr(nrow1)
		mean22=dblarr(nrow2)
		for i=0,nrow1-1 do mean21[i]=mean(d1[i,*])
		for i=0,nrow2-1 do mean22[i]=mean(d2[i,*])
		cov21=correlate(d1,/covariance,/double)
		cov22=correlate(d2,/covariance,/double)
		evals21=eigenql(cov21,/double,eigenvectors=evecs21) > (MACHAR(/double)).xmin
		evals22=eigenql(cov22,/double,eigenvectors=evecs22) > (MACHAR(/double)).xmin
		sevecs21=evecs21
		sevecs22=evecs22
       for i=0,n_elements(evals21)-1 do begin
      	    sevecs21[*,i]=sqrt(evals21[i])*evecs21[*,i]
       endfor
       for i=0,n_elements(evals22)-1 do begin
           	sevecs22[*,i]=sqrt(evals22[i])*evecs22[*,i]
       endfor
           	    	    	    
        detcov21=product(evals21)
        invcov21=la_invert(cov21,/double)
        detcov22=product(evals22)
        invcov22=la_invert(cov22,/double)

        a21=0
       	kfac21=dblarr(ncol1) 
        for i=0,ncol1-1 do begin
        	kfac21[i]=transpose(d1[*,i]-mean21)#invcov21#(d1[*,i]-mean21)
        	a21=a21+kfac21[i]
        endfor
		kmax21=max(kfac21)
		
        a22=0
        kfac22=dblarr(ncol2)
        for i=0,ncol2-1 do begin                        
        	kfac22[i]=transpose(d2[*,i]-mean22)#invcov22#(d2[*,i]-mean22)                 
			a22=a22+kfac22[i]
        endfor
        kmax22=max(kfac22)
   
        ellvol21=ellipsoid_vol(detcov21,max(kfac21),nrow1)
        ellvol22=ellipsoid_vol(detcov22,max(kfac22),nrow2)
        
        ellipsoid1={cov:cov21,invcov:invcov21,detcov:detcov21,mean:mean21,evals:evals21,evecs:evecs21,sevecs:sevecs21,ksum:a21,kmax:kmax21,ellvol:ellvol21}
        ellipsoid2={cov:cov22,invcov:invcov22,detcov:detcov22,mean:mean22,evals:evals22,evecs:evecs22,sevecs:sevecs22,ksum:a22,kmax:kmax22,ellvol:ellvol22}

		;;DETERMINE BAYESIAN INFORMATION CRITERION FOR TWO CLUSTERS. TAKEN FROM MULTINEST BY F. FEROZ
        BIC2=(ncol1+ncol2)*nrow1*alog(2.0*!DPI)+ncol1*alog(detcov21)+ncol2*alog(detcov22)+a21+a22 $
        -2.0*(ncol1*alog(ncol1)+ncol2*alog(ncol2))+2.0*(ncol1+ncol2)*alog(ncol1+ncol2)+(nrow1^2.0 + 3.0*nrow1+1.0)*alog(ncol1+ncol2)
        
    ENDELSE
    
    volfrac=double(ellvol21+ellvol22)/double(ellvol1)
 
	crit1=((BIC2 le BIC1) or (volfrac le 0.0D)) ;; VOLFRAC IS EXPERIMENTAL BUT NOT NECESSARY 
	crit2=(nwh1 gt nrow+1 and nwh2 gt nrow +1)

 	IF crit1 and crit2 THEN BEGIN
 		branch1=cluster_data(d1,ellipsoid=ellipsoid1)
 		branch2=cluster_data(d2,ellipsoid=ellipsoid2)
 		return_val=[branch1,branch2]
 	ENDIF else begin
 		data_pointer=ptr_new(data)
 		return_val={data_ptr:data_pointer,data_num:n_elements(data[0,*]),mean:mean1,cov:cov1,invcov:invcov1,evecs:evecs1,sevecs:sevecs1,evals:evals1,detcov:detcov1,$
 					ellvol:ellipsoid_vol(detcov1,kmax1,nrow),kmax:kmax1,ksum:a}
 	ENDELSE
 	
 	return, return_val
END
