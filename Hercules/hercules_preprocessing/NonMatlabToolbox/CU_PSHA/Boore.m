function PGA = Boore(d,h,m,n)
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
    Coef=[-0.04713 0.6909 0.01130 -0.00202 0.19 0.24 0.29];
    sigma=0.27;
    Coef2=[-3.700012 1.1169 0.00615 -0.00045 0.1 0.25 0.36];
    sigma2=0.30;
elseif n == 1
    Coef=[2.991 0.03525 0.00759 -0.00206 0.019 0.25 0.29];
    sigma=0.23;
    Coef2=[2.301 0.02237 0.00012 0.0 0.1 0.25 0.36];
    sigma2=0.36;
else
    warning('Unrecognized source type')
    return;
end
Delta=0.00724*10^(0.507*m); 
R=(d^2+Delta^2)^0.5;

PGA=[(Coef(1)+Coef(2)*m+Coef(3)*h+Coef(4)*R+Coef(5)*0+Coef(6)*0+Coef(7)*0) sigma...
    (Coef2(1)+Coef2(2)*m+Coef2(3)*h+Coef2(4)*R+Coef2(5)*0)+Coef(6)*0+Coef(7)*0 sigma2];