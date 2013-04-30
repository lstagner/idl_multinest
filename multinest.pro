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
;				   PLOT = plot				)
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

FUNCTION multinest,LOG_LIKELIHOOD_FUNC,PRIOR_FUNC,NUM_PARAMS,NUM=NUM,SAMPLE_NUM=SAMPLE_NUM,TOL=TOL,EXPAND=EXPAND,PLOT=PLOT
	
	DefSysV, '!RNG', Obj_New('RandomNumberGenerator')
	if not keyword_set(NUM) then NUM=1000
	if not keyword_set(SAMPLE_NUM) then SAMPLE_NUM=1
	if not keyword_set(EXPAND) then EXPAND=1.0d
	if not keyword_set(TOL) then TOL=0.001d
	
	;;Get initial samples
	live_points = !RNG -> GetRandomNumbers(NUM_PARAMS,NUM,/double)
	
	live_samples=CALL_FUNCTION(PRIOR_FUNC,live_points)
	print,"GOT INITIAL LIVE SAMPLES"
	;;Get initial likelihoods and likelihood min
	live_likelihoods=dblarr(NUM)
	for i=0L, NUM-1 do begin
		live_likelihoods[i]=CALL_FUNCTION(LOG_LIKELIHOOD_FUNC,live_samples[*,i])
		print,"GOT",i,"th INITIAL LIVE LIKELIHOOD ",live_likelihoods[i]
	endfor
	print,"GOT ALL INITIAL LIVE LIKELIHOODS"	
	likelihood_min=MIN(live_likelihoods,likelihood_min_loc,MAX=likelihood_max)
	sample_likelihoods=live_likelihoods

	;;Initialize Varibles	
	H = 0.0
	logZ=-(machar(/double)).xmax
	k=0.0
	logw=alog(1.0-exp(-1.0/double(NUM)))
	deltaZ = TOL
	
	;;Begin Nested Sampling Loop
	while (deltaZ ge TOL) do begin
	
		;;Cluster the live points 		
		result=cluster_data(live_points)
		
		;;Determine the number of clusters
		clusters=n_elements(result)

		;;Determine the number of points in each cluster
		ellnum=dblarr(clusters)
		for i=0,clusters-1 do begin
			ellnum[i]=result[i].data_num
		endfor
				
		totnum=total(ellnum)
		cnt=0L

		while cnt le SAMPLE_NUM-1 do begin
			;;Pick an ellipse to sample from based on the number of points in the ellipses
			ellprob=total(double(ellnum)/double(totnum),/cumulative)
			r=!RNG->GetRandomDigits(1)
			w=in_ellipsoids(live_points[*,(r mod NUM)],result, expand=EXPAND, /location)
			if w[0] eq -1 then continue
			if n_elements(w) gt 1 then w[0]=w[r mod n_elements(w)]
			cov=result[w[0]].cov
			mean=result[w[0]].mean
			sevecs=result[w[0]].sevecs ;Scaled Eigenvectors                        
			k_fac_max=result[w[0]].kmax ;Ellipsoid surface constant i.e. (x/a)^2 + (y/b)^2 + ... == kmax
			
			;;Sample ellipsoid until a viable point is obtained viable being greater than likelihood min
			while 1 do begin
				samp=sample_ellipsoid(cov,mean,k_fac_max,expand=EXPAND,scaled_evecs=sevecs)
				w=where(samp gt 1.0D or samp lt 0.0D,nw)
				if nw ne 0 then continue
				tsamp=CALL_FUNCTION(PRIOR_FUNC,transpose(samp))
				samp_likelihood=CALL_FUNCTION(LOG_LIKELIHOOD_FUNC,transpose(tsamp))
				if (samp_likelihood gt likelihood_min) then begin
	        		break
				endif
			endwhile
			
			;;Determine which ellipsoids does the sample lie in. Ellipsoids could intersect so even though
			;;we sampled from one ellipsoid it could lie in more than one
			samp_loc=in_ellipsoids(transpose(samp),result,expand=EXPAND)
			r=!RNG->GetRandomNumbers(1,/double)
			;;Accept point with probability 1.0/(number of ellipsoids the point lies in)
			if r le (1.0d/samp_loc) then begin                                                    	    				
				if k eq 0 then samples=tsamp else samples=[[samples],[tsamp]]
				;;Add point to live points
				live_points=[[live_points],[transpose(samp)]]
				live_likelihoods=[live_likelihoods,samp_likelihood]
				
				;;Update Log Evidence(logZ) and Information(H) 
				;;Formula in D.S. Sivia Data Analysis: A Bayesian Tutorial
				logwt=double(logw+likelihood_min)
				logZnew=log_plus(logZ,logwt) 
				H = exp(logwt-logZnew)*likelihood_min + exp(logZ-logZnew)*(H+logZ)-logZnew
				logZ = logZnew
				logw = logw-1.0/(double(NUM))
				
				;;Determine Probability of sample 
				if k eq 0 then log_prob=logwt else log_prob=[log_prob,logwt]
				
				;;Remove the live point with lowest likelihood(likelihood_min)
				w=where(live_likelihoods ne likelihood_min,w_num,complement=nw,ncomplement=nw_num)
				if nw_num gt 1 then w=[w,nw[1:*]]
			    live_likelihoods=live_likelihoods[w]
				live_points=live_points[*,w]
				
				;;Determine newest likelihood min and max
				likelihood_min=MIN(live_likelihoods,likelihood_min_loc,MAX=likelihood_max)
				cnt=cnt+1

			endif
		endwhile
		k=k+cnt
		;;Determine the largest contribution to the log evidence from the existing live points if the sampling loop
		;;If this is less then some tolerance the sampling loop is terminated
		deltaZ = log_plus(logZ,likelihood_max-k/double(NUM))-logZ
		if k mod 1 eq 0 then print,"DELTA LOGZ: ",deltaZ

		if keyword_set(PLOT) then begin
			if NUM_PARAMS eq 1 then !p.multi=0
			if NUM_PARAMS lt 7 then !p.multi=[0,NUM_PARAMS,1,0,0]
			if NUM_PARAMS ge 7 and NUM_PARAMS lt 13 then !p.multi=[0,6,2,0,0]
            for i=0,NUM_PARAMS-1 do begin
                hist=histogram(samples[i,*],locations=xbins,NBINS=50)
                plot,xbins,hist/total(hist),psym=10
            endfor                                                       
		endif
		for i=0,clusters-1 do ptr_free,result[i].data_ptr


	endwhile
	
	;;Add the log evidence contribution from the remaining live points 
	logw=-double(k)/double(NUM)-alog(double(NUM))

	for i=0, NUM-1 do begin
		logwt=double(logw+live_likelihoods[i])		
		logZnew=log_plus(logZ,logwt)
		log_prob=[log_prob,logwt]
		H = exp(logwt-logZnew)*live_likelihoods[i] + exp(logZ-logZnew)*(H+logZ)-logZnew
		logZ=logZnew	
		samples=[[samples],[CALL_FUNCTION(PRIOR_FUNC,live_points[*,i])]]
	endfor
	
	;;Normalize log probability
	log_prob=log_prob-logZ
	;;Sort the points for convienence 
	sort_order=SORT(log_prob)
	samples=samples[*,sort_order]
	log_prob=log_prob[sort_order]
	
	print,"H = " , H 
	print,"Log(Z) = ", logZ , " +- ", sqrt(H/(1.0*NUM))
	
	;;Return data structure that contains the samples,the log probability,
	;;the log evidence and error, and the information respectively
	data_str={samples:samples,log_prob:log_prob,logZ:logZ,logZ_err:sqrt(H/(1.0*NUM)),H:H}
	return, data_str
	GET_OUT:
END
