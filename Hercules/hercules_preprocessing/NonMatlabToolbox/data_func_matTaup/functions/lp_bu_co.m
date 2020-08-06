% lp_bu_co.m
% low pass butterworth filter implementation
% usage:
% filtered_data = lp_bu_co(data,corner (Hz), samples_per_second, order, passes)
% designed to mimic sac
% for example, to low pass filter data at 10 seconds in sac you would type:
% lp bu co 0.1
% with the default order (number of poles) at 2
% this function would have the syntax if the data was at 40 samples per second:
% filtered = lp_bu_co(data,0.1,40,2);
% matlab uses the command butter as:
% [a,b] = butter(order,w) where the frequency is defined as nyquist * w (w is
% from 0 to 1). 
%
% passes should be 1 for causal or 2 for acausal

function y = lp_bu_co(x,corner_f,sps,n,t)
    nyquist = sps / 2;
    w = corner_f / nyquist;
    [b,a] = butter(n,w);
    if t == 1
        y = filter(b,a,x);
    else
        y = filtfilt(b,a,x);
    end
    
    return