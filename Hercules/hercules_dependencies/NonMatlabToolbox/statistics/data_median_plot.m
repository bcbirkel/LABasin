function [dataYInterval medianintervals] = data_median_plot(dataX,dataY,deltaInterval)
%
% DATA_MEAN_PLOT
%

yInterval=min(dataY):deltaInterval:max(dataY);
j=1;
for i=1:length(yInterval)-1    
    a=yInterval(i)
    b=yInterval(i+1)
    k=find(dataY>=a & dataY<=b)
    yN(j)=(a+b)*.5;
    dataXN(j)=median(dataX(k));  
    clear('k','a')
    j=j+1;
end

k=isnan(dataXN);
dataYInterval=yN(~k);
medianintervals=dataXN(~k);

figureR;
hold
plot(dataX,dataY,'k.','LineStyle','none');
plot(medianintervals,dataYInterval);


return
