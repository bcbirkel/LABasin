% hp_bu_co.m
% high pass butterworth filter implementation
% usage:
% filtered_data = hp_bu_co(data,corner (Hz), samples_per_second, order, passes)
% designed to mimic sac
% for example, to high pass filter data at 1 second in sac you would type:
% hp bu co 1
% with the default order (number of poles) at 2
% this function would have the syntax if the data was at 40 samples per second:
% filtered = hp_bu_co(data,1,40,2);
% matlab uses the command butter as:
% [a,b] = butter(order,w,'type') where the frequency is defined as nyquist * w (w is
% from 0 to 1). 
%
% the last argument, passes should be 1 for a causal filter and 2 for an
% acausal filter

function y = hp_bu_co(x,corner_f,sps,n,t)
    nyquist = sps / 2;
    w = corner_f / nyquist;
    [b,a] = butter(n,w,'high');
    if t==1
        y = filter(b,a,x);
    else
        y = filtfilt(b,a,x);
    end
    
    return