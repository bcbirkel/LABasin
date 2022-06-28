%
%   s=surface_corrector_polygon(polygon_x,polygon_y,xgrid,ygrid,zgrid,value)
%
%   input:
%         polgyon   - the polygon that we want to remove and replace with a 
%                     value in a given surface
%         value     - the value that will be put where the polygon is located
%         constrain - if elev is < constrain do not modifiy elevation.
%
%   output: 
%         s       - corrected surface
%
function s=surface_corrector_polygon_from_SurfA_to_surfB(polygon_x,polygon_y,xgrid,ygrid,zgridA,zgridB,constrain)

[x,y]=meshgrid(xgrid,ygrid);
xi=reshape(x,length(xgrid)*length(ygrid),1);
yi=reshape(y,length(xgrid)*length(ygrid),1);
zi=reshape(zgridB,length(xgrid)*length(ygrid),1);
in=inpolygon(xi,yi,polygon_x,polygon_y);
zi(in)=zgridA(in);
k=find(zi<constrain);
zi(k)=zgridB(k);
     
s=reshape(zi,length(ygrid),length(xgrid));

return

