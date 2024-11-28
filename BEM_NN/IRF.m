function[B33] = IRF (ndata,nirf, t, b33, del_t,omega)
B33=zeros(nirf,1);
B33(1) = trapz(omega, b33);
for i=2:nirf
    t(i)=(i-1)*del_t;
    B33(i)=0;
    for j=2:ndata
        sum1 = (b33(j)-b33(j-1))/(omega(j)-omega(j-1));
        sum2 = cos(omega(j)*t(i)) - cos(omega(j-1)*t(i));
        
        B33(i)=B33(i)+ sum1*sum2;
    end
    fact1 = 2.0/(pi()*t(i)^2);
    fact2 = 2.0/(pi()*t(i));
    sum3 = b33(ndata)*sin(omega(ndata)*t(i));
    B33(i)=fact1*B33(i)+ fact2*sum3; %2*pi()*b33(47)*sin(omega(47)*(n-1)*del_t)/t;
end

end
