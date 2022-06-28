%
% density function: Nafe-Drake equation
%
%          input  : vp velocity in km/s         
%          returns: density in kg/m^3
%
function val=density(vp)    
    
    rho= 1.6612*vp-.4721*vp^2+.0671*vp^3-.0043*vp^4+.000106*vp^5

    val=rho*1000;

return

