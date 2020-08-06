% coherence.m
% compute the magnitude squared coherence:
% Cfg * Cfg / (Cff * Cgg)
% from two time series (f(t) and g(t))
%
% usage:
% out = coherence(f,g,spsf,spsg)
%
% where f and g are the timeseries's and spsf and spsg are the samples per
% second for f and g respectively.
% decimates if necessary so the two functions have the same sampling rate.
%
% also, checks the lengths of f and g and assumes the same start time by
% only keeping the number of points present in the shorter time period
%
% matlabs dsp internal program mscohere is similar, but it instead applies
% a welch average to compute psd's which does not return exactly what is
% always desired. This is designed to be more simple..

%DO NOT USE

function [out] = coherence(f,g,spsf,spsg)
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
    end
    
    % convert to frequency domain
    ff=fft(f);
    gg=fft(g);
    
    % compute the Gxx, Gxy, and Gyy and put into a temporary array
%    numer = ((real(ff) .* real(gg) + imag(ff) .* imag(gg)) + (imag(ff) .* real(gg) + real(ff) .* imag(gg)));
%    denom = complex_amplitude(ff) .* complex_amplitude(gg);
%    tmp = numer ./ denom;
%   t1 here is actually wrong. A welch's average is necessary or else this will always properly compute to 1.
	t1 = real(ff) .* real(ff) + imag(ff) .* imag(ff);


	t2 = real(gg) .* real(gg) + imag(gg) .* imag(gg);
    t3 = real(ff) .* real(gg) + imag(gg) .* imag(ff);
    tmp = (t3 .* t3) / (t1 .* t2);
	
    % kill half spectra
    out = tmp(1:length(tmp)/2);
    
    % cap the coherence s it goes from 0 to 1
    for j=1:length(out)
        if out(j) > 1
            out(j) = 1;
        end
        if (out(j) < 0)
            out(j) = 0;
        end
    end
    
    return
