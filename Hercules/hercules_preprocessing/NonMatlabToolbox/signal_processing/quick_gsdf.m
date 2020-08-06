function GSDF = quick_gsdf(Syn,Obs,Iso,dt,Frec)
%
% GSDF = gsdf(Syn,Obs,IsoF,dt,Freq)
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
%   GSDF - Matrix (m*5) containig the gsdf measurements [F dTp dTq]
%

% GSDF filtering parameters
Tw = 2/mean(Frec);      % Effective window
Sw = 2*pi*0.72/Tw;      % Sigma w
NT = length(Syn);       % Number of samples
LagT = (-dt*(NT-1):dt:dt*(NT-1))';      % Lag time
T = 0:dt:dt*(NT-1);

% Lag Time can't be larger than 5% length of the signal
W = tapercos(LagT,-0.1*T(end),0.1*T(end),0.2);

% normalize to avoid small numbers
Obs = Obs/max(Obs);
Syn = Syn/max(Syn);
Iso = Iso/max(Iso);

% % % figure;
% % % subplot(3,2,1),
% % % plot(T,Obs,'linewidth',2), hold on
% % % plot(T,Syn,'linewidth',1.5)
% % % plot(T,Iso,'linewidth',1.5)
% % % ylim(1.1*[-1 1])

% Crosscorrelations of the Signals
ObsC(:,1) = xcorr(Iso,Obs,'biased');
SynC(:,1) = xcorr(Iso,Syn,'biased');

% Find max to center Window
[~,InO] = max(W.*ObsC);
[~,InS] = max(W.*SynC);

% Windowing
WO = exp(-(0.5*Sw^2)*(LagT-LagT(InO)).^2).*ObsC;
WS = exp(-(0.5*Sw^2)*(LagT-LagT(InS)).^2).*SynC;  

% % % subplot(3,2,2),
% % % plot(LagT,ObsC/max(SynC),'linewidth',2), hold on
% % % plot(LagT,SynC/max(SynC),'linewidth',1.5)
% % % plot(LagT,exp(-(0.5*Sw^2)*LagT.^2),'linewidth',1.5)
% % % ylim([-2 2])
% % % 
% % % subplot(3,2,3), 
% % % plot(LagT,WO/max(WS),LagT,WS/max(WS),'linewidth',2), hold on
% % % plot(LagT,exp(-(0.5*Sw^2)*LagT.^2),'linewidth',1.5)
% % % ylim([-2 2])


% Parameters for "bootstraping"
[~,InOR] = max(WO); % Maximo de las correlaciones
[~,InSR] = max(WS);

% Allocate outputs
GSDF = zeros([length(Frec) 3]);

% Compute for all frequencies
for nF = 1:length(Frec)
    
    Si=0.1*2*pi*Frec(nF);   % Sigma i
    
    % Crosscorrelagram and Autocorrelagram filtering
    dO = computebandfftfilter_gauss(WO',dt,Frec(nF),Si);
    dS = computebandfftfilter_gauss(WS',dt,Frec(nF),Si);
    
% % %     subplot(3,2,4)
% % %     plot(LagT,dO/max(dO),'linewidth',2), hold on
% % %     plot(LagT,dS/max(dO),'linewidth',1.5), hold on
    
    % Let's do bootstraping
    [~,InO]=max(dO); 
    [~,InS]=max(dS);
    
    BS = 1; Cn = 0;
    while BS == 1 || Cn < 10
        if (LagT(InO) < LagT(InOR)+0.51/Frec(nF)) && (LagT(InO) >= LagT(InOR)-0.51/Frec(nF))
            BS = 0;
        elseif LagT(InO) >= LagT(InOR)+0.45/Frec(nF)
                InO=InO-round(1/Frec(nF)/dt);
        elseif LagT(InO) < LagT(InOR)-0.45/Frec(nF)
                InO=InO+round(1/Frec(nF)/dt);
        end
        Cn = Cn+1;
    end
    
    BS = 1; Cn = 0;
    while BS == 1 || Cn < 10
        if (LagT(InS) < LagT(InSR)+0.5/Frec(nF)) && (LagT(InS) >= LagT(InSR)-0.5/Frec(nF))
            BS = 0;
        elseif LagT(InS) >= LagT(InSR)+0.45/Frec(nF)
                InS=InS-round(1/Frec(nF)/dt);
        elseif LagT(InS) < LagT(InSR)-0.45/Frec(nF)
                InS=InS+round(1/Frec(nF)/dt);
        end
        Cn = Cn+1;
    end

    % Save the GSDF measurements
    GSDF(nF,1) = Frec(nF);
    GSDF(nF,2) = LagT(InO)-LagT(InS);	% delta_P
    GSDF(nF,3) = -log(dO(InO)/dS(InS))/(2*pi*Frec(nF));  % Amplitud
    
end