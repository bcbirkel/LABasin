function [GSDF,Jpq,PS,PO]=gsdf_simple(Syn,Obs,IsoF,dt,Frec)
%
% [GSDF,Jp,PO,PS] = gsdf(Syn,Obs,IsoF,dt,Freq)
% Computes de Generalized Seismological Data Functionals following Gee &
% Jordan (1992)
%
% Inputs:
%   Syn - Synthetic seismogram (n*1)
%   Obs - Observed waveform (n*1)
%   IsoF - Isolation filter (n*1)
%   dt - Time step of all the prev signals signals
%   Frec - Vector with the frequencies for the measurements (m*1)
%
% Outputs:
%   GSDF - Matrix (m*5) containig the gsdf measurements [F dTp dTg dTq dTq]
%   Jpq - Complex perturbation kernel for dTp and dTq (Cheng et al. (2007))
%   PS and PO - Vectors with the five parameters of Gaussian wavelet
%               approximation of the Filtered Windowed Crosscorrelograms
%
% This code uses the nlinfit  function for fitting the five parameter 
% Gaussian wavelet, check results if warnings are shown.

% GSDF filtering parameters
Tw = 2/mean(Frec);       % Effective window
sw = 2*pi*0.72/Tw;      % Sigma w
no = length(Syn);       % Number of samples
tt = (-dt*(no-1):dt:dt*(no-1))';  % Lag time
T = 0:dt:dt*(no-1);

% figure;
% subplot(3,2,1),
% plot(T,Obs/max(Obs),'linewidth',2), hold on
% plot(T,1.5+Syn/max(Syn),'linewidth',1.5)
% plot(T,3+IsoF/max(Syn),'linewidth',1.5)
% xlim([3 7])

% Crosscorrelations of the Signals
ObsC(:,1)=xcorr(IsoF,Obs,'biased');
SynC(:,1)=xcorr(IsoF,Syn,'biased');

% A local maximum will be selected if 0.25<CC<4 and closest to 0-lag
[PkS,InS] = findpeaks(envelope(SynC,1,'peak'),'MinPeakHeight',0.25);
[PkO,InO] = findpeaks(envelope(ObsC,1,'peak'),'MinPeakHeight',0.25);

if isempty(InS) || isempty(InO)
    [~,InS] = max(SynC);
    [~,InO] = max(ObsC);
else
    [~,In] = min((1+abs(tt(InS))).^2.*abs(1-PkS));
    InS = InS(In);

    [~,In] = min((1+abs(tt(InO))).^2.*abs(1-PkO));
    InO = InO(In);
end

% Windowing
WO = exp(-(0.5*sw^2)*(tt-tt(InO)).^2).*ObsC;
WS = exp(-(0.5*sw^2)*(tt-tt(InS)).^2).*SynC;  


% subplot(3,2,2), plot(tt,ObsC/max(SynC),'linewidth',2), hold on
% plot(tt,SynC/max(SynC),'linewidth',1.5)
% plot(tt,exp(-(0.5*sw^2)*(tt-tt(InS)).^2),'linewidth',1.5)
% xlim([-2 2])
% 
% subplot(3,2,3), plot(tt,WO/max(WS),tt,WS/max(WS),'linewidth',2), hold on
% plot(tt,exp(-(0.5*sw^2)*(tt-tt(InS)).^2),'linewidth',1.5)
% xlim([-2 2])


% Parameters for "bootstraping"
[~,InOR] = max(WO); % Maximo de las correlaciones
[~,InSR] = max(WS);

% Model for fitting
Eqn = @(b,x)(b(1)*exp(-0.5*b(2).^2*(x-b(3)).^2).*cos(b(4)*(x-b(5))));
ne = min(round(2/min(Frec)/dt),round(no/2));    % Effective bandwidth for inversion

%Parametros de GSDF
w0=2*pi/(tt(end));               % Frecuencias
wN=2*pi/(2*dt);
w(:,1)=-wN:w0:wN;

% Isolation filter FFT for perturbation kernel
tff=conj(fftshift(fft(IsoF)))*1/no;
tff = reshape(tff,no,1);

% Allocate outputs
GSDF = zeros([length(Frec) 5]);
PS = zeros([length(Frec) 5]);
PO = zeros([length(Frec) 5]);
Jpq = zeros([length(Frec) no]);

% Compute for all frequencies
for nF = 1:length(Frec)
    
    si=0.1*2*pi*Frec(nF);   % Sigma i
    % Crosscorrelagram and Autocorrelagram filtering
    dO=computebandfftfilter_gauss(WO',dt,Frec(nF),si);
    dS=computebandfftfilter_gauss(WS',dt,Frec(nF),si);
    
%     subplot(3,2,4)
%     plot(tt,1.5*nF+dO/max(dO),'linewidth',2), hold on
%     plot(tt,1.5*nF+dS/max(dO),'linewidth',1.5), hold on
%     xlim([-2 2])
    
    % Let's do bootstraping
    [~,InO]=max(dO); 
    [~,InS]=max(dS);
    
    BS = 1; Cn = 0;
    while BS == 1 || Cn < 10
        if (tt(InO) < tt(InOR)+0.51/Frec(nF)) && (tt(InO) >= tt(InOR)-0.51/Frec(nF))
            BS = 0;
        elseif tt(InO) >= tt(InOR)+0.45/Frec(nF)
                InO=InO-round(1/Frec(nF)/dt);
        elseif tt(InO) < tt(InOR)-0.45/Frec(nF)
                InO=InO+round(1/Frec(nF)/dt);
        end
        Cn = Cn+1;
    end
    
    BS = 1; Cn = 0;
    while BS == 1 || Cn < 10
        if (tt(InS) < tt(InSR)+0.5/Frec(nF)) && (tt(InS) >= tt(InSR)-0.5/Frec(nF))
            BS = 0;
        elseif tt(InS) >= tt(InSR)+0.45/Frec(nF)
                InS=InS-round(1/Frec(nF)/dt);
        elseif tt(InS) < tt(InSR)-0.45/Frec(nF)
                InS=InS+round(1/Frec(nF)/dt);
        end
        Cn = Cn+1;
    end

    % Five parameter Gaussian wavelet fitting
    [As,Is] = max(envelope(dS));
    [Ao,Io] = max(envelope(dO));
    GaS = nlinfit(tt(Is-ne:Is+ne),dS(Is-ne:Is+ne),Eqn,[As 2*pi*si tt(Is) 2*pi*Frec(nF) tt(InS)]);
    GaO = nlinfit(tt(Io-ne:Io+ne),dO(Io-ne:Io+ne),Eqn,[Ao 2*pi*si tt(Io) 2*pi*Frec(nF) tt(InO)]);
    
    % Check fitting
    if ((GaO(1)/GaS(1)) > 1e5) || abs(GaO(5)-GaS(5)) > tt(end)/2
        GaO = [Ao 2*pi*si tt(Io) 2*pi*Frec(nF) tt(InO)];
        GaS = [As 2*pi*si tt(Is) 2*pi*Frec(nF) tt(InS)];
    end
     
%     subplot(2,1,1)
%     plot(tt,Eqn(GaS,tt),tt,dS);
%     subplot(2,1,2)
%     plot(tt,Eqn(GaO,tt),tt,dO);
%     close;

    % Parametros:
    wi=2*pi*Frec(nF);
    
    wP=((si^2)*w+(sw^2)*wi)/(sw^2+si^2);
    wPP=((si^2)*w-(sw^2)*wi)/(sw^2+si^2);
    siP=((si^2)*(sw^2)/(sw^2+si^2))^0.5;
    
    % Estimate seismogram perturbation kernel
    IW=(siP/(sw*GaS(1)))*exp(-0.5*(w-2*pi*Frec(nF)).^2/(sw^2+si^2)).*tff./wP+...
        (siP/(sw*GaS(1)))*exp(-0.5*(w+2*pi*Frec(nF)).^2/(sw^2+si^2)).*tff./wPP;
        
    IW(1:floor(size(IW)/2))=0*IW(1:floor(size(IW)/2));
    
    itff = ifft(fftshift(no*IW));
    
    % Save the GSDF measurements
    GSDF(nF,1) = Frec(nF);
    GSDF(nF,2) = GaO(5)-GaS(5);	% delta_P
    GSDF(nF,3) = GaO(3)-GaS(3);	% delta_G
    GSDF(nF,4) = -log(GaO(1)/GaS(1))/GaS(4);  % Amplitud
    GSDF(nF,5) = -(GaO(4)-GaS(4))/GaS(2)^2;  % Frequency
    
    %Save the seismogram perturbation kernel
    itff(isnan(itff)) = 0;
    Jpq(nF,:) = itff;
    
    % Save the Gaussian approximation parameters
    PO(nF,:) = GaO;
    PS(nF,:) = GaS;
    clear GaO GaS InO InS;
end