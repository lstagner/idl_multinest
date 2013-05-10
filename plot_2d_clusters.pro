PRO plot_2d_clusters,clusters,x_range=x_range,y_range=y_range
 	if not keyword_set(x_range) then x_range=[0,1]
 	if not keyword_set(y_range) then y_range=[0,1]
    num_clusters=n_elements(clusters)

    plot,[0,0],[1,1],xrange=x_range,yrange=y_range,/nodata
    for i=0,num_clusters-1 do begin
       	loadct,13,/silent
       	mean=clusters[i].mean
       	cov=clusters[i].cov
       	invcov=clusters[i].invcov
       	evals=clusters[i].evals
       	evecs=clusters[i].evecs

       	d_ptr=clusters[i].data_ptr
       	d=*d_ptr
       	k=dblarr(clusters[i].data_num)

       	for j=0L, clusters[i].data_num-1 do begin
            tmp=d[*,j]-mean
            k[j]=tmp##invcov##transpose(tmp)
        endfor

        kmax=max(k,nk)
        for j=0,n_elements(evals)-1 do begin
            evecs[*,j]=sqrt(evals[j])*evecs[*,j]
        endfor

        rd=dblarr(2,1000)
        for j=0L, n_elements(rd[0,*])-1 do begin
            rd[*,j] = sample_ellipsoid(cov,mean,kmax,expand=1.1,scaled_evecs=evecs,/surface)
        endfor

        oplot,d[0,*],d[1,*],psym=3,color=50+i*50
		loadct,0,/silent
        oplot,rd[0,*],rd[1,*],psym=3

    endfor

END
