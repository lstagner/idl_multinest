FUNCTION kmedoids, data,n_clusters=n_clusters,n_iterations=n_iterations

	if not keyword_set(n_clusters) then n_clusters=2L
	if not keyword_set(n_iterations) then n_iterations=10
	
	n_points=n_elements(data[0,*])
	distances=DISTANCE_MEASURE(data,/double,/matrix)
	w=where(distances ne 0)
	energies=distances
;	energies[w]=1.0/distances[w]
	if n_clusters eq 2L then begin
		tmp=max(distances,w)
		ind=array_indices(distances,w)
		centers=data[*,ind]
	endif else begin
		r = !RNG->GetRandomNumbers(n_clusters,/long)
		ind= r mod n_points
		centers=data[*,ind]
	endelse

	energy=energies[ind,*]
	tmp=min(energy,tmp2,dimension=1)
	tmp3=array_indices(energy,tmp2)
	cluster_id=transpose(tmp3[0,*])
	tot_energy=total(energy[cluster_id,lindgen(n_points)])

	for i=0,n_iterations-1 do begin
		r=!RNG->GetRandomNumbers(1,/long)
        k = r mod n_points
		tot_energy_new=dblarr(n_clusters)
		for j=0,n_clusters-1 do begin
			ind_new=ind
			if k eq ind_new[j] then continue
			ind_new[j]=k
			energy_new=energies[ind_new,*]
    		tmp=min(energy_new,tmp2,dimension=1)
			tmp2=array_indices(energy_new,tmp2)
			cluster_id_new=transpose(tmp2[0,*])
			tot_energy_new[j]=total(energy_new[cluster_id_new,lindgen(n_points)])
		endfor
		w=where(tot_energy_new lt tot_energy,nw)
		if nw ne 0 then begin
			tmp=min(tot_energy_new,w)
			ind[w]=k
			centers=data[*,ind]
		endif
	endfor

	return, centers
END
