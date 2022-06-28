%
%    brochervelfit: Brocher's regression fit
%
%          input  : vs velocity in km/s         
%          returns: vp in km/s
%
function val=brochervelfitvpgivenvs(vs)    
    
    val = 0.9409+2.0947*vs-0.8206*vs^2+0.2683*vs^3-.0251*vs^4; 
    

return
