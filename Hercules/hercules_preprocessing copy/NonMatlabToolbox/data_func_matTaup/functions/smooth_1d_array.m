% smooth_1d_array
% usage:
% ss = smooth_1d_array(input_array,nsmooth)
%
% smoothes the array, input_array, by nsmooth points and returns the array
% into ss
function ss = smooth_1d_array(in,n)
    % initialize output for speed
    ss=zeros(1,length(in));
    
    % check cases whether odd or even
    if mod(n,2) == 1
        % smooth over the input array
        for j=1:length(in)
           % initialize for the averaging
           sum=0;
           count=0;
           % check the extremes so we don't run out of bounds
           if j-floor(n/2) < 1
               start = 1;
           else
               start = j-floor(n/2);
           end
           if j+floor(n/2)>length(in)
                ee = length(in);
           else
               ee = j+floor(n/2);
           end
           % now do the smoothing 
           for jj=start:ee
               sum = sum + in(jj);
               count = count + 1;
           end
           ss(j) = sum / count;
        end
    else
        for j=1:length(in)
           % initialize for the averaging
           sum=0;
           count=0;
           % check the extremes so we don't run out of bounds
           if j-floor(n/2) < 1
               start = 1;
           else
               start = j-floor(n/2);
           end
           if j+floor(n/2)-1>length(in)-1
                ee = length(in)-1;
           else
               ee = j+floor(n/2)-1;
           end
           % now do the smoothing 
           for jj=start:ee
               sum = sum + (in(jj)+in(jj+1))/2;
               count = count + 1;
           end
           ss(j) = sum / count;
        end
        
    end
    
    return
