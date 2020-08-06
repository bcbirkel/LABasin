function PGA = gmpeArr(d,m)
% Strong Ground-Motion Attenuation Relationship (PGA)
% PGA = gmpeRP(d,h,m) 
% Input:
%   d - hypocentral distance
%   m - magnitud
% Output:
%   PGA - 1x2 vector with PGA and sigma
% Reference: Arroyo et al., 2010.

Coef=[2.4862 0.9392 0.5061 0.0150 -0.3850 -0.0181];
Coef2=[-7.1389 1.8721 0.5376 0.0001 0.5592 -0.05334];

r0=(1.4447e-5*exp(2.3026*m))^0.5;
RR=(d^2+r0^2)^0.5;

E1=expint(Coef(4)*d);

E2=expint(Coef(4)*RR);

E12=expint(Coef2(4)*d);

E22=expint(Coef2(4)*RR);

PGA=[Coef(1)+Coef(2)*m+Coef(3)*log((E1-E2)/r0^2) 0.75 ...
    Coef2(1)+Coef2(2)*m+Coef2(3)*log((E12-E22)/r0^2) 0.6701];