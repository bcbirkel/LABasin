%
%  Help:
%
% [f absfourier,phasefourier,realfourier,imagfourier,fourier]=computefft(signal,dt)
%
%

function [f absfourier,phasefourier,realfourier,imagfourier,fourier]=computefft(signal,dt)


timeSteps = size(signal);

if(timeSteps(1)==1)
    numsteps=timeSteps(2);
end

if(timeSteps(2)==1)
    numsteps=timeSteps(1);
end

 t=0:dt:(numsteps-1)*dt;


fourier=fft(signal);

f=(1/dt)*( 0:(numsteps-1) )/numsteps;
absfourier=abs(fourier);
phasefourier=angle(fourier);
realfourier=real(fourier);
imagfourier=imag(fourier);
fourier=dt*fourier;
return;
