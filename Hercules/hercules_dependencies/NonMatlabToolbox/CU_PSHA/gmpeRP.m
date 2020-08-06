function PGA = gmpeRP(d,h,m,n)
% Strong Ground-Motion Attenuation Relationship (PGA)
% PGA = gmpeRP(d,h,m) 
% Input:
%   d - hypocentral distance
%   h - hypocentral depth
%   m - magnitud
%   n - near trench n=1 or inslab n=0
% Output:
%   PGA - 1x2 vector with PGA and sigma
% Reference: Rodriguez-Perez, 2014.

if n == 0
    Coef=[-1.1321 0.8038 0.0033 -0.0014 1.3219];
    sigma=0.39;
    Coef2=[-3.7248 0.9135 0.0015 -0.0014 0.8602];
    sigma2=0.32;
elseif n == 1
    Coef=[-1.2324 0.5016 0.0141 -0.0006 0.9432];
    sigma=0.037;
    Coef2=[-2.7903 0.5844 -0.0050 0.0008 0.7850];
    sigma2=0.34;
else
    warning('Unrecognized source type')
    return;
end
Delta=0.00724*10^(0.507*m); 
R=(d^2+Delta^2)^0.5;

PGA=[(Coef(1)+Coef(2)*m+Coef(3)*h+Coef(4)*R-Coef(5)*log10(R)) sigma...
    (Coef2(1)+Coef2(2)*m+Coef2(3)*h+Coef2(4)*R-Coef2(5)*log10(R)) sigma2];