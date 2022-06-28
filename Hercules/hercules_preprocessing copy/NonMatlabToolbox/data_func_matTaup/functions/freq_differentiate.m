% freq_differentiate.m
% performs differentiation in the frequency domain by applying an fft,
% multiplying by 2*pi*f, and then ifft back to time domain
%
% usage:
% out = freq_differentiate(data,samples_per_second)
% designed to differentiate from displacement to velocity or velocity to
% acceleration
%
% see companion freq_integrate to integrate from acceleration to
% velocity or velocity to displacement

% advantage of this method is that frequency domain is more stable than
% time domain and thus you don't get large offsets. It may be suggested to
% do a high pass filter after just to remove long period offset due to the
% lack of integration constant

function out = freq_differentiate(data,sps)
    a=fft(data);
    b=zeros(1,length(data));
    for j=1:length(data)
        omega = 2*pi*(sps/2/length(data)) * j;
        b(j) = real(a(j)) * omega * i - imag(a(j)) * omega;
    end
    out = ifft(b,'symmetric');
    return