%
%    brochervelfit: Brocher's regression fit
%
%          input  : vp velocity in km/s         
%          returns: vs in km/s
%
function val=brochervelfitvsgivenvp(vp)    
    
    val = 0.7858-1.2344*vp+0.7949*vp^2-0.1238*vp^3+.0064*vp^4 
    

return
