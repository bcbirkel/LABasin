%
%    brochervelfit: Brocher's regression fit
%
%          input  : vp velocity in km/s         
%          returns: vs in km/s
%
function [vpbr,vsbr,rhobr]=computebrochersproperties(vp,vs,defaultvp)    

if(vp == 0 && vs==0 )
    vpbr  =defaultvp;                    
    rhobr =density(vpbr); 
    vsbr  =brochervelfitvsgivenvp(vpbr);  
elseif(vs == 0)
    vpbr  =vp;
    rhobr =density(vp);
    vsbr  =brochervelfitvsgivenvp(vpbr);
elseif(vp == 0)
    vpbr  =brochervelfitvpgivenvs(vs);                    
    rhobr =density(vpbr);
    vsbr  =vs;
else
    vpbr = vp;
    vsbr = vs;
    rhobr= density(vpbr);
end 

vpbr=vpbr*1000;
vsbr=vsbr*1000;

% Check poisson ratio
% Check lambda value
% Check mu value 


return
