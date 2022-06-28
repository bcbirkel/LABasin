%
%   levels_struct=surface_extract_levels_given_Z(xgrid,ygrid,zgrid,level_array)
%
%   input:
%
%   output: 
%
function levels_struct=surface_extract_levels_given_Z(xgrid,ygrid,zgrid,levels,smooth)

if(smooth==1)
    zgrid=surface_filter(zgrid,2,2)
end    

c=contour(xgrid,ygrid,zgrid,levels)

level(1)=c(1,1);
numPoints(1)=c(2,1);
index(1)=2;
index(2)=numPoints(1)+2;
% obtain the different index
i=2;
indexValid=1;

while (indexValid == 1)
    level(i)=c(1,index(i));
    numPoints(i)=c(2,index(i));
    index(i+1)=c(2,index(i))+index(i)+1;
    contourToPlot(i-1).values=c(:,index(i)+1:index(i+1)-1);
        
    if(index(i+1) > length(c))
        indexValid=-1
    end
   i=i+1;
    
end


%figureR;
%hold
% Plot the contours and make a stuct
for i=1:length(numPoints)-1
    line(contourToPlot(i).values(1,:),contourToPlot(i).values(2,:) ,contourToPlot(i).values(2,:)*0+level(i))
    levels_struct(i).level=level(i);
    levels_struct(i).x=contourToPlot(i).values(1,:);
    levels_struct(i).y=contourToPlot(i).values(2,:);
end
