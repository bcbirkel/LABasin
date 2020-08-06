%
%   surfaceF=surface_filter_constant(surface,filterSize,type)
%
%   input:
%               type 1-average
%                    2-gaussian
%
%   output:
%
function surfaceF=surface_filter(surface,filterSize,type)

if(type==1)
    h=ones(filterSize,filterSize)/(filterSize*filterSize);
    s=size(surface);
    surfaceF=filter2(h,surface);
else
    
    n1=filterSize;sigma1=filterSize;n2=filterSize/2;sigma2=filterSize/2;theta1=0;
    filt=d2gauss(n1,sigma1,n2,sigma2,theta1);
    filt=filt/sum(sum(filt));
    surfaceF=filter2(filt,surface);
           
end

for i=1:filterSize
    surfaceF(i,:)=surfaceF(filterSize,:);
    surfaceF(:,i)=surfaceF(:,filterSize);
end


return

