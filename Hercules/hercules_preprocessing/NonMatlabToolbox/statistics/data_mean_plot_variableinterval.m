function [dataYInterval meanintervals stdinterval] = data_mean_plot(dataX,dataY,yInterval)
%
% DATA_MEAN_PLOT
%

%yInterval=min(dataY):deltaInterval:max(dataY);
j=1;
for i=1:length(yInterval)-1    
    a=yInterval(i)
    b=yInterval(i+1)
    k=find(dataY>=a & dataY<=b)
    yN(j)=(a+b)*.5;
    % remove nan
    dataAll=dataX(k);
    clear k
    k=isnan(dataAll)
    dataXN(j)=mean(dataAll(~k));  
    dataXNSD(j)=std(dataAll(~k));
    clear('k','a','dataAll')
    j=j+1;
end

k=isnan(dataXN);
dataYInterval=yN(~k);
meanintervals=dataXN(~k);
stdinterval=dataXNSD(~k);

% figureR;
% hold
% plot(dataX,dataY,'k.','LineStyle','none');
% plot(meanintervals,dataYInterval);
% 
% 
return
