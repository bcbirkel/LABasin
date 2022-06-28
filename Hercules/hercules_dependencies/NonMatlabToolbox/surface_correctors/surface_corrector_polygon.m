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
function s=surface_corrector_polygon(polygon_x,polygon_y,xgrid,ygrid,zgrid,value)

[x,y]=meshgrid(xgrid,ygrid);
xi=reshape(x,length(xgrid)*length(ygrid),1);
yi=reshape(y,length(xgrid)*length(ygrid),1);
zi=reshape(zgrid,length(xgrid)*length(ygrid),1);
%in=inpolygon(xi,yi,polygon_x,polygon_y);
%check for nans
nanindex=isnan(polygon_x);

% do it polygon by polygon (if they are nested)
separators=[1 find(nanindex==1)]
if(separators==1)
    separators(2)=length(polygon_x)+1
end
for i=1:length(separators)-1
    if(i==1)
        indexA=1;
    else
        indexA=separators(i)+1
    end
    indexB=separators(i+1)-1;
    polx=polygon_x(indexA:indexB)'
    poly=polygon_y(indexA:indexB)'
    
    dimensions=size(poly)
    if(dimensions(1)==1)
        polx=polx';
        poly=poly';
    end
    
    in=inpoly([xi yi],[polx poly]);
    zi(in)=value;   
end

s=reshape(zi,length(ygrid),length(xgrid));

return

