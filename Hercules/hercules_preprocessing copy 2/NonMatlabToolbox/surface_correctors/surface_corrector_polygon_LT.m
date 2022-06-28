%
%   s=surface_corrector_polygon(polygon_x,polygon_y,xgrid,ygrid,zgrid,value)
%
%   input:
%         polgyon - the polygon that we want to remove and replace with a 
%                   value in a given surface
%         value   - the value that will be put where the polygon is located
%
%   output: 
%         s       - corrected surface
%
function s=surface_corrector_polygon_LT(polygon_x,polygon_y,xgrid,ygrid,zgrid,value)

[x,y]=meshgrid(xgrid,ygrid);
xi=reshape(x,length(xgrid)*length(ygrid),1);
yi=reshape(y,length(xgrid)*length(ygrid),1);
zi=reshape(zgrid,length(xgrid)*length(ygrid),1);
in=inpolygon(xi,yi,polygon_x,polygon_y);

ziN=zi;
ziN(in)=value;
zi=max(zi,ziN);

s=reshape(zi,length(ygrid),length(xgrid));

return

