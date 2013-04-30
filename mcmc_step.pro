PRO mcmc_step,start,constraint,log_likelihood_func,prior_func,samp,err,NUM_STEPS=NUM_STEPS

;	DefSysV, '!RNG', Obj_New('RandomNumberGenerator')
	;;SET NUMBER OF MCMC STEPS
	if not keyword_set(NUM_STEPS) then NUM_STEPS=10
	COMMON mcmc_common

	tstart = CALL_FUNCTION(prior_func,start)
	p = CALL_FUNCTION(log_likelihood_func,transpose(tstart))
	x = start
	r = !RNG -> GetRandomNumbers(n_elements(start),NUM_STEPS,/double,/NORMAL)
	i=0
	while i lt NUM_STEPS do begin
		tot+=1.0d
		x_prime = x + sigma*r[*,i]
;		print,x_prime,x
		w=where(x_prime gt 1.0D or x_prime lt 0.0D,nw)

		if nw ne 0 then continue
		xt_prime=CALL_FUNCTION(prior_func,transpose(x_prime))
		p_prime=CALL_FUNCTION(log_likelihood_func,transpose(xt_prime))
		a = p_prime/p
		if (r[i] le a) and (p_prime gt constraint) then begin
			x = x_prime
			p = p_prime
			accept+=1.0d
		endif else tot+=1.0d
		i = i+1 
	endwhile	

	if accept/tot gt .9 then sigma = sigma * exp(1.0d / accept)
	if accept/tot lt .6 then sigma = sigma * exp(-1.0d / (tot-accept))

	samp=x
	if array_equal(transpose(samp),start) eq 1 then err=1
	samp=transpose(samp)
END 