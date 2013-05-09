;NAME:
;	SAMPLE_ELLIPSOID
;AUTHOR:
;	LUKE STAGNER
;	UNIVERSITY OF CALIFORNIA, IRVINE
;	LSTAGNER@UCI.EDU
;PURPOSE:
;	UNIFORMALLY SAMPLES WITHIN AN ELLIPSOID GIVEN THE MEAN, C0VARIANCE MATRIX AND K FACTOR OF THE ELLIPSOID
;CALLING SEQUENCE:
;	POINT = SAMPLE_ELLIPSOID(	COV,						$
;								MEAN,						$
;								KFAC,						$
;								EXPAND=EXPAND,				$
;								SCALED_EVECS = SCALED_EVECS,$
;								SURFACE = SURFACE			)
;
;PARAMETERS:
;	COV: COVARIANCE MATRIX OF THE ELLIPSOID 
;	MEAN: MEAN OF THE ELLIPOID
;	KFAC: K FACTOR OF ELLIPOID I.E. (X/A)^2 + (Y/B)^2 + ... == KFAC 
;KEYWORDS:
;	EXPAND: EXPANSION FACTOR OF K FACTOR
;	SCALED_EVECS: EIGENVECTORS OF THE COVARIANCE MATRIX SCALED BY SQRT OF THE RESPECTIVE EIGENVALUES
;	SURFACE: IF SET THE RETURN VALUE WILL BE A RANDOM POINT ON THE SURFACE OF THE ELLIPSOID
;RETURN VALUE:
;	A N-DIMENSIONAL ARRAY THAT CONTAINS A RANDOMLY SAMPLED POINT WITHIN THE ELLISOID
;REQUIRED PROGRAMS:
;	RANDOMNUMBERGENERATOR__DEFINE.PRO
;COMMON BLOCKS:
;	NONE
;REFERENCES: 
;	MULTINEST PROGRAM BY F. FEROZ
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION sample_ellipsoid,cov,mean,kfac,expand=expand,scaled_evecs=scaled_evecs,surface=surface,ptr=ptr

	;;IN SCALED_EVECS EIGENVECTORS ARE IN ROWS
;    DefSysV, '!RNG', Obj_New('RandomNumberGenerator')
    
    ;;SET INITIAL VALUES
	if not keyword_set(expand) then expand=1.0D
	if not keyword_set(NUM) then NUM=1	
 
 	dim=double(n_elements(mean))	
	
	;;DETERMINE SCALED EIGENVECTORS
	if not keyword_set(scaled_evecs) then begin
		evals=EIGENQL(cov,/double,eigenvectors=evecs)
		scaled_evecs=dblarr(n_elements(evals),n_elements(evals))
		for i=0,n_elements(evals)-1 do begin
			scaled_evecs[*,i]=sqrt(evals[i])*evecs[*,i]
		endfor 
	endif	
	
	;;GET NORMALLY DISTRIBUTED RANDOM SAMPLE
	rn=!RNG->GetRandomNumbers(fix(dim),1,/NORMAL)
	rn=double(rn)
	if not keyword_set(surface) then begin 
		ru=!RNG->GetRandomNumbers(1)
		ru=double(ru[0])
	endif else ru=1.0D
	
	;;TRANSFORM NORMALLY SAMPLED POINT TO A UNIFORM POINT IN THE ELLIPOIDS
	;;TAKEN FROM MULTINEST BY F. FEROZ
	tmp= (ru^( 1.0D / dim )) / sqrt( total( rn^2.0D ) ) 
	rn=tmp[0]*rn
	samp= sqrt(expand*kfac) * transpose(scaled_evecs)##transpose(double(rn)) + mean
	if keyword_set(ptr) then return,ptr_new(reform(samp)) else return,samp
END
