% power_spectrum.m
% function to compute a power spectrum
% To compute a full power spectral density, use this along with average and
% smoothing code to extract a reasonable mean daily psd estimate
% usage:
% spect = power_spectrum(ts,sps)
% For best results, preprocess the data by removing a mean, linear trend,
% applying a cosine taper, and removing the instrument response

function aa = power_spectrum(ts,dt)
    % fft the input array
    bb = fft(ts);
    
    % kill half spectra
    cc=bb(1:length(ts)/2);
    
    % power = 10 * log10((real*real+imag*imag)*2*dt)/(npts)
    aa=10*log10((real(cc).*real(cc)+imag(cc).*imag(cc))*dt*2/length(ts));
    
    
    return
