
function IFFT = RedNoise(N)
% IFFT = RedNoise(N) creates a time serieswhich spectrum is corrsponds to
% red noise. The time series has standard deviation around 1

    nhaf = floor(N/2);
    freq = [0:nhaf]/(N-1);

    lag1 = 0.72;

    Real = (1-lag1^2) ./ (1-2*lag1*cos(freq*2*pi)+lag1^2);  % [Eqn(16)]

    Phi = 2*pi*rand([1 nhaf+1]);
    FFT = Real.*exp(1j*Phi);


    IFFT = 10*real(ifft([0 FFT(1:end) flip(FFT(2:end-1))]));
end