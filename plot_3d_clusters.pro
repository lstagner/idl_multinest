PRO plot_3d_clusters,clusters
    num_clusters=n_elements(clusters)

	LOADCT, 39, NCOLORS=!D.N_COLORS-1,/silent
	TVLCT, 255, 255, 255, !D.N_COLORS-1
	; Set the 3D coordinate space with axes.
	SURFACE, DIST(5), /NODATA, /SAVE, XRANGE=[0.0,1.0], $
	         YRANGE=[0.0,1.0], ZRANGE=[0.0, 1.0], XSTYLE=1, $
	         YSTYLE=1, ZSTYLE=1, CHARSIZE=1.5, COLOR=255, $
    	     BACKGROUND=0 ;;change az,ax to change viewing angle
       
	for i=0,num_clusters-1 do begin
	    d_ptr=clusters[i].data_ptr
        d=*d_ptr
        d=d[0:2,*]

		PLOTS,d[0,*],d[1,*],d[2,*], PSYM=2, COLOR=50+i*50, /T3D
		PLOTS,d[0,*],d[1,*],d[2,*]*0,psym=3,color=50+i*50,/T3D
		PLOTS,d[0,*],d[1,*]*0+1,d[2,*],psym=3,color=50+i*50,/T3D
		PLOTS,d[0,*]*0+1,d[1,*],d[2,*],psym=3,color=50+i*50,/T3D
		
	endfor
END
