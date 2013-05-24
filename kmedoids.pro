FUNCTION kmedoids, data,n_clusters=n_clusters,n_iterations=n_iterations

	if not keyword_set(n_clusters) then n_clusters=2L
	if not keyword_set(n_iterations) then n_iterations=10
	
	n_points=n_elements(data[0,*])
	distances=DISTANCE_MEASURE(data,/double,/matrix)

	if n_clusters eq 2L then begin
		tmp=max(distances,w)
		ind=array_indices(distances,w)
		centers=data[*,ind]
	endif else begin
		r = !RNG->GetRandomNumbers(n_clusters,/long)
		ind= r mod n_points
		centers=data[*,ind]
	endelse

	cost=distances[ind,*]
	tmp=min(cost,tmp2,dimension=1)
	tmp2=array_indices(cost,tmp2)
	cluster=transpose(tmp2[0,*])
	tot_cost=total(cost[cluster,lindgen(n_points)])		

	for i=0,n_iterations-1 do begin
		r=!RNG->GetRandomNumbers(1,/long)
        k = r mod n_points
		tot_costp=dblarr(n_clusters)
		for j=0,n_clusters-1 do begin
			indp=ind
			if k eq indp[j] then continue
			indp[j]=k
			costp=distances[indp,*]
    		tmp=min(cost,tmp2,dimension=1)
			tmp2=array_indices(cost,tmp2)
			cp=transpose(tmp2[0,*])
			tot_costp[j]=total(cost[cluster,lindgen(n_points)])
		endfor
		w=where(tot_costp lt tot_cost,nw)
		if nw ne 0 then begin
			tmp=min(tot_costp,w)
			ind[w]=k
			centers=data[*,ind]
		endif
	endfor

	return, centers
END
