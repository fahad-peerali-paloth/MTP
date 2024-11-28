function [ Cradij ] = radiation_restoration(Aijw,k,w,T,Aijinf)
n = length(T);

if iscolumn(k) == 0
    k = k';
end
if iscolumn(T) == 0
    T = T';
end

if n>100
    k(101:n) = k(100)*ones(n-100,1);
end
% a*k +b
if n ~=1
    a = (k(2:n) - k(1:n-1))./(T(2:n) - T(1:n-1));
    b = ((k(1:n-1).*T(2:n)) - (k(2:n).*T(1:n-1))) ./ (T(2:n) - T(1:n-1));

    int = ((-a/w) .* (T(2:n).*cos(w*T(2:n)) - T(1:n-1).*cos(w*T(1:n-1))));
    int = int + ((a/w^2) .* (sin(w*T(2:n)) - sin(w*T(1:n-1))));
    int = int - ((b/w).*(cos(w*T(2:n)) - cos(w*T(1:n-1))));

    Cradij = w^2*(Aijinf - Aijw) - w*sum(int);
else
    Cradij = 0;
end
end

