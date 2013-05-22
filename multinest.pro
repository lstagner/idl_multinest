;NAME:
;	MULTINEST
;AUTHOR:
;	LUKE STAGNER 
;	UNIVERSITY OF CALIFORNIA, IRVINE
;	LSTAGNER@UCI.EDU
;VERSION: 1.0 
;PURPOSE:
;	MULTINEST PERFORMS MULTIMODEL NESTED SAMPLING
;DESCRIPTION:
;	MULTINEST USES X-MEANS CLUSTERING TO PARTITION THE SET OF LIVE POINTS. EACH CLUSTER IS THEN SURROUNDED BY AN 
;	ELLIPSOID. AN ELLIPSOID IS SELECTED ACCORDING TO THE NUMBER OF LIVE POINTS ENCLOSED. A SAMPLE IS THEN DRAWN 
;	FROM THE ELLIPSOID SUBJECT TO THE LIKELIHOOD CONSTRAINT. THIS IS REPEATED UNTIL A CERTAIN TOLERANCE .  
;CALLING SEQUENCE:
;	S = MULTINEST( LOG_LIKELIHOOD_FUNC,		$
;				   PRIOR_FUNC,				$
;				   NUM_PARAMS,				$
;				   NUM = num,				$
;				   SAMPLE_NUM = sample_num, $
;				   TOL = tol				$
;				   EXPAND = expand,         $
;				   PLOT = plot, 			$	
;				   SILENT=silent			)
;PARAMETERS:
;	LOG_LIKELIHOOD_FUNC: NAME OF LOG LIKELIHOOD FUNCTION OF THE FORM: VALUE=LOG_LIKELIHOOD_FUNC(ARRAY[NUM_PARAMS]) 
;						 (LOG LIKELIHOOD FUNCTION RETURNS THE LOG OF THE LIKELIHOOD FUNCTION GIVEN A SET OF PARAMETERS)
;	PRIOR_FUNC: NAME OF PRIOR FUNCTION OF THE FORM: ARRAY[NUM_PARAMS,NUM] =  PRIOR_FUNC(ARRAY[NUM_PARAMS,NUM])
;				(PRIOR FUNCTION MAPS UNIFORM RANDOM SAMPLES FROM HYPERCUBE TO WANTED PRIOR PROBABILITY DISTRIBUTIONS)
;	NUM_PARAMS: THE NUMBER OF PARAMETERS BEING ESTIMATED
;KEYWORDS:
;	NUM: NUMBER OF LIVE POINTS DEFAULT SET TO 1000
;	SAMPLE_NUM: NUMBER OF SAMPLES DRAWN BETWEEN CLUSTERING DEFAULT SET TO 1
;	TOL: LOG EVIDENCE TOLERANCE DEFAULT SET TO .001 IN LOG EVIDENCE 
;	EXPAND: EXPANSION FACTOR FOR SAMPLING ELLIPSOIDS 
;	PLOT: SHOWS HISTOGRAMS FOR EACH PARAMETER BEING ESTIMATED
;RETURN VALUE:
;	A STRUCTURE WHICH CONTAINS THE SAMPLES DRAWN ACCORDING TO THE POSTERIOR DISTRIBUTION, THE INFORMATION, THE LOG_EVIDENCE, 
;	ASSOCIATED LOG_EVIDENCE ERROR, AND LOG_PROBABILITY FOR EACH SAMPLE
;REQUIRED PROGRAMS:
;	CLUSTER_DATA.PRO
;	SAMPLE_ELLIPSOID.PRO
;	IN_ELLIPSOIDS.PRO
;	RANDOMNUMBERGENERATOR__DEFINE.PRO
;	PLOT_2D_CLUSTERS.PRO
;COMMON BLOCKS:
;	NONE 
;REFERENCES:
;F. FEROZ, M.P. HOBSON, "MULTIMODAL NESTED SAMPLING: AN EFFICIENT AND ROBUST ALTERVATIVE TO MCMC METHODS FOR ASTRONOMICAL DATA ANALYSIS"
;	Mon. Not. Roy. Astron. Soc., 384, 2, 449-463 (2008)
;
;F. FEROZ, M.P. HOBSON, M. BRIDGES, "MULTINEST: AN EFFICIENT AND ROBUST BAYESIAN INFERENCE TOOL FOR COSMOLOGY AND PARTICLE PHYSICS" 
;	Mon. Not. Roy. Astron. Soc. 398: 1601-1614 (2009)
;DATE MODIFIED:
;	08/06/2012

FUNCTION multinest,LOG_LIKELIHOOD_FUNC,PRIOR_FUNC,NUM_PARAMS,$
				   NUM=NUM,SAMPLE_NUM=SAMPLE_NUM,TOL=TOL,EXPAND=EXPAND,$
				   PLOT=PLOT,SILENT=SILENT,NTHREADS=NTHREADS,USE_THREADS=USE_THREADS
	
	DefSysV, '!RNG', Obj_New('RandomNumberGenerator')
	if not keyword_set(NUM) then NUM=1000
	if not keyword_set(SAMPLE_NUM) then SAMPLE_NUM=1
	if not keyword_set(EXPAND) then EXPAND=1.0d
	if not keyword_set(TOL) then TOL=0.001d
	if not keyword_set(NTHREADS) then NTHREADS=4
	
	if keyword_set(USE_THREADS) then begin
		;;Create threads
		create_threads,NTHREADS,threads,/idl_startup
		execute_threads,threads,"DefSysV, '!RNG', Obj_New('RandomNumberGenerator') & !EXCEPT=0"
		wait_threads,threads
	endif

	COMMON mcmc_common,accept,tot,sigma
	accept=0.0d & tot=0.0d  & sigma=.01d

	;;Get initial samples
	live_points=ptrarr(NUM,/allocate_heap)
	for i=0,NUM-1 do begin
		*live_points[i] = !RNG -> GetRandomNumbers(NUM_PARAMS,/double)
	endfor

	live_samples=CALL_FUNCTION(PRIOR_FUNC,live_points)

;	print,"GOT INITIAL LIVE SAMPLES"
	;;Get initial likelihoods and likelihood min
	live_likelihoods=dblarr(NUM)
	for i=0L, NUM-1 do begin
		live_likelihoods[i]=CALL_FUNCTION(LOG_LIKELIHOOD_FUNC,live_samples[i])
	endfor

	likelihood_min=MIN(live_likelihoods,likelihood_min_loc,MAX=likelihood_max)
	
	;;Initialize Varibles	
	H = 0.0
	logZ=-(machar(/double)).xmax
	k=0L
	logw=alog(1.0-exp(-1.0/double(NUM)))
	deltaZ = TOL

	samples=ptrarr(10000)
	log_prob=dblarr(10000,/nozero)
	n_samp=10000L	
	;;Begin Nested Sampling Loop
	while (deltaZ ge TOL) do begin
;		print,'Clustering Data'
		;;Cluster the live points 		
		result=cluster_data(live_points)
		
		;;Determine the number of clusters
		clusters=n_elements(result)
		print,'Number of clusters:',clusters
		;;Determine the number of points in each cluster
		ellnum=dblarr(clusters)
		for i=0,clusters-1 do begin
			ellnum[i]=result[i].data_num
		endfor
		totnum=total(ellnum)
		cnt=0L

		if keyword_set(plot) then begin
			if NUM_PARAMS ne 2 then no_ellipse=1
			plot_2d_clusters,result,no_ellipse=no_ellipse
		endif
		index=0
		while cnt le SAMPLE_NUM-1 do begin
			hit=0.0 & miss=0.0
			;;Pick an ellipse to sample from based on the number of points in the ellipses
			ellprob=total(double(ellnum)/double(totnum),/cumulative)
			r=!RNG->GetRandomNumbers(1,/long)
			start_point=live_points[(r mod NUM)]
			w=in_ellipsoids(live_points[(r mod NUM)],result, expand=EXPAND, /location)
			if w[0] eq -1 then continue
			r=!RNG->GetRandomNumbers(1,/long) ;;pick new random number to be independent of the choice of ellipsoid
			if n_elements(w) gt 1 then w[0]=w[r mod n_elements(w)]
			cov=result[w[0]].cov
			mean=result[w[0]].mean
			sevecs=result[w[0]].sevecs ;Scaled Eigenvectors                        
			k_fac_max=result[w[0]].kmax ;Ellipsoid surface constant i.e. (x/a)^2 + (y/b)^2 + ... == kmax
			if keyword_set(USE_THREADS) then begin
				in_vars=['PRIOR_FUNC','LOG_LIKELIHOOD_FUNC','cov','mean','k_fac_max','EXPAND','sevecs']
				out_vars=['samp','tsamp','samp_likelihood']
				push_variable,threads,in_vars
				wait_threads,threads
			endif
			;;Sample in Parallel
			samp_likelihood=likelihood_min-1
			count=0L
			while 1 do begin
				if count lt 100 then begin
					if keyword_set(USE_THREADS) then begin
						if index eq 0 then begin
							cmd='get_sample,PRIOR_FUNC,LOG_LIKELIHOOD_FUNC,cov,mean,k_fac_max,EXPAND,sevecs,samp,tsamp,samp_likelihood'
							execute_threads,threads,cmd
							wait_threads,threads
							pull_variable,threads,out_vars,output
						endif
						tmp=index
						for i=tmp,NTHREADS-1 do begin
							if clusters eq 1 then index+=1
							samp_likelihood= (*(output[i,2]))
							if samp_likelihood gt likelihood_min then begin
								samp=output[i,0]
								tsamp=output[i,1]
								if clusters eq 1 then print,index
								break
							endif
						endfor
						if index ge NTHREADS then index=0
					endif else begin
						get_sample,PRIOR_FUNC,LOG_LIKELIHOOD_FUNC,cov,mean,k_fac_max,EXPAND,sevecs,samp,tsamp,samp_likelihood,/ptr
					endelse
					count+=1
				endif else begin
;					print,'mcmc start'
					mcmc_step,start_point,likelihood_min,LOG_LIKELIHOOD_FUNC,PRIOR_FUNC,samp,tsamp,samp_likelihood,err,/ptr
					if err eq 1 then begin
						samp_likelihood=likelihood_min-1.
;						print,'mcmc failed'
					endif; else print,'mcmc finished'
					count=0L
				endelse
				if samp_likelihood gt likelihood_min then break else miss=miss+1
			endwhile		

			;;Determine which ellipsoids does the sample lie in. Ellipsoids could intersect so even though
			;;we sampled from one ellipsoid it could lie in more than one
			samp_num=in_ellipsoids(samp,result,expand=EXPAND)
			r=!RNG->GetRandomNumbers(1,/double)

			;;Accept point with probability 1.0/(number of ellipsoids the point lies in)
			if r le (1.0d/samp_num) then begin
				hit+=1.0

				;;Check if we need more space                            	    				
				if (k+cnt) ge n_samp then begin
					 samples=[samples,ptrarr(NUM)]
					 log_prob=[log_prob,dblarr(NUM,/nozero)]
					 n_samp+=NUM	 
				endif
				
				;;Determine newest likelihood min and max
				new_likelihood_min=MIN([live_likelihoods,samp_likelihood])
	
				;;Update Log Evidence(logZ) and Information(H) 
				;;Formula in D.S. Sivia Data Analysis: A Bayesian Tutorial
;				print,'Updating Evidence'
				logwt=double(logw + log_plus(new_likelihood_min,likelihood_min) - alog(2.0d));;Trapazoid Integration
				logZnew=log_plus(logZ,logwt) 
				H = exp(logwt-logZnew)*likelihood_min + exp(logZ-logZnew)*(H+logZ)-logZnew
				logZ = logZnew
				logw = logw-1.0d/(double(NUM))

				;;add to samples and record prob 
				log_prob[k+cnt]=logwt
				samples[k+cnt]=live_samples[likelihood_min_loc]

				;;Add point to live points
				live_points[likelihood_min_loc]=samp
				live_samples[likelihood_min_loc]=tsamp
				live_likelihoods[likelihood_min_loc]=samp_likelihood
				
			    likelihood_min=MIN(live_likelihoods,likelihood_min_loc,MAX=likelihood_max)

				cnt=cnt+1
			endif else miss+=1
		endwhile
		k=k+cnt
		;;Determine the largest contribution to the log evidence from the existing live points if the sampling loop
		;;If this is less then some tolerance the sampling loop is terminated
		deltaZ = log_plus(logZ,likelihood_max-(k+cnt)/double(NUM))-logZ
		if k mod SAMPLE_NUM eq 0 and not keyword_set(SILENT) then print,'DELTA LOGZ: '+strtrim(string(deltaZ))+string(hit/(hit+miss))

		for i=0,clusters-1 do ptr_free,result[i].data_ptr
	endwhile
	
	;;Add the log evidence contribution from the remaining live points 
	logw=-double(k)/double(NUM)-alog(double(NUM))

	for i=0, NUM-1 do begin
		logwt=double(logw+live_likelihoods[i])		
		logZnew=log_plus(logZ,logwt)
		H = exp(logwt-logZnew)*live_likelihoods[i] + exp(logZ-logZnew)*(H+logZ)-logZnew
		logZ=logZnew
    	if (k+cnt) ge n_samp then begin
        	samples=[samples,ptrarr(NUM)]
			log_prob=[log_prob,dblarr(NUM,/nozero)]
			n_samp+=NUM
        endif
		log_prob[k]=logwt
        samples[k]=live_samples[i]
		k+=1
	endfor
	samples=samples[0:k-1]
	log_prob=log_prob[0:k-1]
	;;Normalize log probability
	log_prob=log_prob-logZ
	;;Sort the points for convienence 
	sort_order=SORT(log_prob)
	samples=samples[sort_order]
	log_prob=log_prob[sort_order]
	
	print,'H = ' , H 
	print,'Log(Z) = ' + string(logZ) + ' +- '+ string(sqrt(H/(1.0*NUM)))
	
	;;Return data structure that contains the samples,the log probability,
	;;the log evidence and error, and the information respectively
	data_str={samples:samples,log_prob:log_prob,logZ:logZ,logZ_err:sqrt(H/(1.0*NUM)),H:H}
	if keyword_set(USE_THREADS) then destroy_threads,threads
	return, data_str
	GET_OUT:
END
