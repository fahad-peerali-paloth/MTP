function[damp] = DampingCoeff(i,nirf, B33, v, del_t)
damp=0;
%t(i) = (i-1)*del_t;
    %prediction
%     z(i) = z(i-1)+v(i-1)*del_t;
%     clear temp_trapz
    if(i<nirf)
        k=i;
    else
        k=nirf;
    end
    
    for j = 1:k-1
        
         damp =  damp + del_t* B33(j)*v(i-j);
     
    end


end