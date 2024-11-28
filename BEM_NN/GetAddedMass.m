function[add] = GetAddedMass(omega, del_t,b33,a33,ndata)  
   
    add = 0;
    
    
    sum1 = 1.0/(omega(ndata)*omega(ndata));
    sum2 = (1.0/(omega(ndata)))*(b33(1) - b33(ndata)*cos(omega(ndata)*ndata*del_t));
    sum3=0;
    for i = 2 : ndata
        t = (i-1)*del_t;
        sum3 = sum3 + ((b33(i) - b33(i-1))/del_t)*(sin(omega(ndata)*t) - sin(omega(ndata)*(t-del_t)));

    end
    
    sum = sum1*sum3 + sum2;
        
    add = a33(ndata) + (1/omega(ndata))*sum;
    
    
end

    



