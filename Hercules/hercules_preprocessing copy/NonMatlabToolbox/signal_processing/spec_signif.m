% [signif] = spec_signif(d,dt,nr)
%
% Calcula el nivel de significancia espectral de la senal d (d es un vector
% [1 x n]) con intervalos de muestreo dt o de una senal con varianza d, con
% respecto a un nivel de ruido rojo de nr [0-1]

function [signif,fft_theor] = spec_signif(Y,siglvl)


N = length(Y);
nhaf = floor(N/2);

freq = [0:nhaf]/(N-1);

variance = var(Y);

lag1 = 0.72;

% get the appropriate parameters [see Table(2)]
fft_theor = (1-lag1^2) ./ (1-2*lag1*cos(freq*2*pi)+lag1^2);  % [Eqn(16)]
fft_theor = siglvl*variance*fft_theor;  % include time-series variance
signif = fft_theor;

return


