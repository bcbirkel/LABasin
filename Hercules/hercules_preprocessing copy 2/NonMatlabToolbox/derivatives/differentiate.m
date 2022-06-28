%
% differentiate(signal)
%
%                 using a three points
%                 need to divide by dt to get the correct value scaled
%

function y=differentiate(signal)

    % k=1 & m=2 See Abramowitz and Stegun
    n=length(signal);
    
    factor=1/2;
    
    coeff(1).val=[-3  4 -1]; % sample 1
    coeff(2).val=[-1  0  1]; % any other
    coeff(3).val=[ 1 -4  3]; % sample n

    y=[];
    %sample i and sample n
    y(1)=dot(coeff(1).val,[signal(1  ) signal(2  ) signal(3)]);
    y(n)=dot(coeff(3).val,[signal(n-2) signal(n-1) signal(n)]); 
    % any other
    for i=2:n-1
       y(i)=dot(coeff(2).val,[signal(i-1) signal(i) signal(i+1)]);
    end
            
    y=y*factor;
    return

end