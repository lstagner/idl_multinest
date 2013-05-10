PRO get_sample,PRIOR_FUNC,LOG_LIKELIHOOD_FUNC,cov,mean,k_fac_max,expand,sevecs,point,samp,log_likelihood,ptr=ptr
	while 1 do begin
    	pnt=sample_ellipsoid(cov,mean,k_fac_max,expand=EXPAND,scaled_evecs=sevecs,/ptr)
		w=where(*pnt gt 1.0D or *pnt lt 0.0D,nw)
		if nw ne 0 then begin
            continue
        endif else break
	endwhile
    tsamp=CALL_FUNCTION(PRIOR_FUNC,pnt)
    log_likelihood=CALL_FUNCTION(LOG_LIKELIHOOD_FUNC,tsamp)
	if keyword_set(ptr) then begin
		point=pnt
		samp=tsamp[0]
	endif else begin
		point=*pnt
		samp=*tsamp[0]
	endelse
END

	
