FUNCTION log_plus , a , b                  
	if a gt b then begin           
    	tmp = double(a) + alog(1.0+exp(double(b)-double(a)))
    endif else begin        
        tmp = double(b) + alog(1.0+exp(double(a)-double(b)))
    endelse
                            
 	return, tmp
END