% cross_correlate.m
% compute the cross correlation
% from two time series (f(t) and g(t))
%
% usage:
% cc = cross_correlate(f,g,spsf,spsg,lag_time)
%
% Optionally, append 'norm' to the usage to normalize the correlation.
%
%
% where f and g are the timeseries's and spsf and spsg are the samples per
% second for f and g respectively.
% decimates if necessary so the two functions have the same sampling rate.
%
% checks the lengths of f and g and assumes the same start time by
% only keeping the number of points present in the shorter time period
%
% after checking the given information, does an fft of each data stream,
% computes the cross correlation as:
% cc(a+bi,c+di) = (ac+bd) + (bc - ad)i
% does a symmetric ifft and only returns the values from
% [-lag_time,lag_time]


function [cc] = cross_correlate(f,g,spsf,spsg,lag,varargin)
    % check sample rates. Compute a common sample rate and decimate
    if (spsf ~= spsg)
        common_sps = gcd(spsf,spsg);
        sf = spsf/common_sps;
        sg = spsg/common_sps;
        f=decimate(f,sf);
        g=decimate(g,sg);
    end
    
    % check record lengths. Print a warning if they don't match
    if length(f) ~= length(g)
        disp('Warning: the length of the two data files is inconsistent. Assuming common start time and keeping the shorter length only!');
        if length(f) < length(g)
            common_length = length(f);
        else
            common_length = length(g);
        end
        ff=zeros(1,common_length);
        gg=zeros(1,common_length);
        for j=1:common_length
            ff(j)=f(j);
            gg(j)=g(j);
        end
        clear f g
        f=ff;
        g=gg;
        clear ff gg
    else
        common_length = length(f);
    end
    
    % convert to frequency domain
    ff=fft(f);
    gg=fft(g);
    
    % free up some memory
    clear f g
     
    % compute the cross correlation
    tmp = (real(ff) .* real(gg) + imag(ff) .* imag(gg)) + (imag(ff) .* real(gg) - real(ff) .* imag(gg)) * i;
%    tmp = rand(1,common_length);

    % check and normalize if asked
    if nargin == 6
        if strcmp(varargin{1},'norm')
            denom = complex_amplitude(ff) .* complex_amplitude(gg);
            tmp = tmp ./ denom;
        end
    end
    
    % free memory
    clear ff gg
    
    % return to the time domain
    tt = ifft(tmp,'symmetric');
  
    % rearrange the time series
    cc=zeros(1,2*lag+1);
    cc(lag+1) = tt(1);
    for j=2:(lag+1)
        if (common_length - j > 0 && lag-j > 0) 
            cc(lag+j) = tt(j);
            cc(lag-j) = tt(common_length-j);
        end
    end
    
    clear tt tmp
    
    return
