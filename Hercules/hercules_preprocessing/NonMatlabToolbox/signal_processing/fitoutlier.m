function y=fitoutlier(x,t)
% Y=fitoutlier(X,T) This code will just remove the outliers and linearly 
% interpolate over their positions using the closest values that are not 
% outliers.
% X is the signal and t is the threshold. Tipicaly is considered t=3 
% If you assume that outliers are more than three standard deviations from 
% the median, for example.
y=x;
all_idx = 1:length(x);
outlier_idx = abs(x - median(x)) > t*std(x);  % Find outlier idx
y(outlier_idx) = interp1(all_idx(~outlier_idx), y(~outlier_idx), all_idx(outlier_idx)); % Do the same thing for y