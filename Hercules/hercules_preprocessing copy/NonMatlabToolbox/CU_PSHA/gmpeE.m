function PGA = gmpeE(d,m)
% Strong Ground-Motion Attenuation Relationship (PGA)
% PGA = gmpeE(d,h,m) 
% Input:
%   d - hypocentral distance
%   m - magnitud
% Output:
%   PGA - 1x2 vector with PGA and sigma
% Reference: Esteva et al., 1973.

Coef=[5600 0.8 -2.0];
RR=d+40;
PGA=[Coef(1)*exp(Coef(2)*m)*(RR^Coef(3)) 0.64];