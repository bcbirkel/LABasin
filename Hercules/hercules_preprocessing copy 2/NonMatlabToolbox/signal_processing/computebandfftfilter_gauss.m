%% ************************************************************************ 
%
% Compute a Butterworth band-pass filter in frequency
%
% [stfinal]=computebandfftfilter_gauss(signal,dt,fc,sigma)
%
%         signal - signal to filter
%             dt - time step
%            fc - central freq
%            sigma - gaussian sigma
%  
%    It will add the same amount of 0s as the size of the signal before
%    filtering.
%    The lower part of the filter will have 2*p
%
%% ************************************************************************

function [stfinal]=computebandfftfilter_gauss(signal,dt,fc,sigma)
signal=[signal(1)*linspace(1,1,100) signal signal(end)*linspace(1,1,100)]';
fc=fc*2*pi;
sigma=sigma*2*pi;

timeSteps = size(signal);
if(timeSteps(1)==1)
    numsteps=timeSteps(2);
end
if(timeSteps(2)==1)
    numsteps=timeSteps(1);
end
NFFT = 2^nextpow2(length(signal));
fourier=fft(signal,NFFT);

%dt=.5;
f   = 2*pi*(1/dt)*( 0:(NFFT-1) )/NFFT;
fil = (exp(-(f-fc).^2./(2*sigma^2))); 
signfiltf = fil'.*fourier;

% complete the congugate
for ij=1:NFFT/2-1
   signfiltf(NFFT-ij+1)=conj(signfiltf(ij+1));     
end
    
sft = ifft(signfiltf);
%signal=interp1(0:length(sft)-1,sft,time);
stfinal=sft(101:length(signal)-100);
return;
