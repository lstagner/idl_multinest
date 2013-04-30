;NAME:
;	IN_ELLIPSOIDS
;AUTHOR:
;	LUKE STAGNER
;	UNIVERSITY OF CALIFORNIA,IRVINE
;	LSTAGNER@UCI.EDU
;PURPOSE:
;	DETERMINES THE NUMBER ELLIPSOIDS A POINT LIES IN. MADE TO WORK WITH CLUSTER_DATA.PRO
;DESCRIPTION:
;
;CALLING SEQUENCE:
;	ARR = IN_ELLIPSOIDS( POINT,				 $
;						 ELLIPSOIDS,		 $
;						 EXPAND=EXPAND 		 $
;						 LOCATION=LOCATION   )
;
;PARAMETERS:
;	POINT: POINT IN N-DIMENSIONAL SPACE
;	ELLIPOIDS: A STRUCTURE THAT CONTAINS INFORMATION AB0UT THE ELLIPSOIDS. SEE THE RETURN VALUE OF CLUSTER_DATA.PRO
;KEYWORDS:
;	EXPAND: THE EXPANSION FACTOR FOR THE ELLIPSOIDS
;	LOCATION: IF SET RETURNS AN ARRAY WITH THE LOCATIONS OF THE ELLISOIDS THAT THE POINT LIES IN
;RETURN VALUE:
;	RETURNS THE NUMBER OF ELLIPSOIDS THE POINT LIES IN
;REQUIRED PROGRAMS:
;	THE ELLIPSOIDS PARAMETER IS USUALLY AQUIRED BY CLUSTER_DATA.PRO
;COMMON BLOCKS:
;	NONE
;REFERENCES:
;	NONE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION in_ellipsoids,point, ellipsoids,expand=expand,location=location ;;return value of cluster_data

	num_ellipsoids=n_elements(ellipsoids)
	
	if not keyword_set(expand) then expand=1.0D
	num=dblarr(num_ellipsoids)
	for i=0, num_ellipsoids-1 do begin
		mean=double(ellipsoids[i].mean)
		invcov=double(ellipsoids[i].invcov)
		kmax=double(ellipsoids[i].kmax)
		kfac=(point-mean)##invcov##(point-mean)
		if kfac[0] le kmax*expand then num[i]=1.0D
	endfor
	loc=where(num ne 0,nloc)
	if keyword_set(location) then $
		return, loc else $
		return, total(num)
END
	