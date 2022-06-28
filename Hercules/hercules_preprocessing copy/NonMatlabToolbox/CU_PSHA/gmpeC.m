function PGA = gmpeC(d,h,m)
% Strong Ground-Motion Attenuation Relationship (PGA)
% PGA = gmpeC(d,h,m) 
% Input:
%   d - hypocentral distance
%   h- hypocentral depth
%   m - magnitud
% Output:
%   PGA - 1x2 vector with PGA and sigma
% Reference: Crouse et al., 1991.

PGA=[10^(6.36+1.76*m-2.73*log(d+1.58*exp(0.0608*m))+0.00916*h) 0.773];