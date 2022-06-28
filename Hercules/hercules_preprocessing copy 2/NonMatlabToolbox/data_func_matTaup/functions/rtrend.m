% rtrend.m
% function to remove a linear trend like in sac
% simply a wrapper for the matlab function 'detrend'
% usage:
% detrended = rtrend(data);

function y=rtrend(x,o,bp)
    if nargin < 2
    o  = 'linear';
end
if nargin < 3
    bp = 1;
end
isrowx = isrow(x);
if isrowx
    x = x(:);   % If a row, turn into column vector
end


switch o
    case {0,'c','constant'}
        y = bsxfun(@minus,x,mean(x,1));   % Remove just mean from each column
        
    case {1,'l','linear'}
        N = size(x,1);
        bp = unique([1; double(bp(:)); N]);   % Include both endpoints
        bp = bp(bp >= 1 & bp <=N);   % Should error in the future
        lbp = length(bp);
        % Build regressor with linear pieces + DC
        a = zeros(N,lbp,class(x));
        a(1:N,1) = (1:N)./N;
        for k = 2:(lbp-1)
            M = N - bp(k);
            a((bp(k)+1):end,k) = (1:M)./M;
        end
        a(1:N,end) = 1;       
        y = x - a*mldivide(a,x);   % Remove best fit
        
    otherwise
        % This should eventually become an error.
        warning(message('MATLAB:detrend:InvalidTrendType', num2str( o )));
        y = detrend(x,'linear',bp);
end

if isrowx
    y = y.';
end
    
    