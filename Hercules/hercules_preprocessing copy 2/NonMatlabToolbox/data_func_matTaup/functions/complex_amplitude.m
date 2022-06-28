% complex_amplitude.m
% usage:
% amp=complex_amplitude(a+bi)
% simple code to compute the amplitude of a complex value
function amp = complex_amplitude(a)
    amp = sqrt(real(a).*real(a)+imag(a).*imag(a));
    return
    