% y = denoise(x,t) denoises the signal x using the S-transform. The
% denosising method follows Parolai (2009).
% inputs:   x - signal and noise
%           t - noise samples (noise variance is estimated from the signal
%           samples 1-t)
% output:   y - denoised signal

function y = denoise(s,t)

ns = length(s);
s = reshape(s,1,ns);

% S-transform noise
tn = stran([s(1:t) 0*s(t+1:end)],1);

% Frequency dependent threshold parameters
Sr = std(real(tn),0,2);	Sr = repmat(Sr,[1 ns]);

Lr = Sr*(2*log(ns))^0.5;    
Gr = 0.5*Lr;
Al = 0.7;

% S transform of the signal
ts = stran(s,1);
xr = real(ts);
xi = imag(ts);


C1 = (abs(xr) >= Lr);
C2 = (abs(xr) <= Gr);

Fxr = (Al*Lr.*((abs(xr)-Gr)./(Lr-Gr)).^2).*((Al-3)*((abs(xr)-Gr)./(Lr-Gr))+4-Al);

Fxr(C1) = xr(C1)-(1-Al)*sign(xr(C1)).*Lr(C1);
Fxr(C2) = 0;

C3 = (abs(xi) >= Lr);
C4 = (abs(xi) <= Gr);

Fxi = (Al*Lr.*((abs(xi)-Gr)./(Lr-Gr)).^2).*((Al-3)*((abs(xr)-Gr)./(Lr-Gr))+4-Al);

Fxi(C3) = xi(C3)-(1-Al)*sign(xi(C3)).*Lr(C3);
Fxi(C4) = 0;


y = istransform(Fxr+1j*Fxi);
