function [wf,sf] = compute_wf_sf_nlf(Syn,dt,Freq)
% [wf,sf]=compute_wf_sf(Syn,,dt,Frec)
% Computes de wf and sf parameters according to Gee & Jordan (1992).
% Inputs are the Synthetic waveform Syn, dt is the time step of the signal
% and Frec is a vector with the frequencies for the fitting;

% Crosscorrelation and autocorrelaion
no=length(Syn);
tt(:,1)=-dt*(no-1):dt:dt*(no-1);

SynC(:,1)=xcorr(Syn,Syn,'biased');

% Pick the maximum of the correlation
Tw=2/mean(Freq);
sw=2*pi*0.72/Tw;    % Sigma w

[~,InS] = max(SynC);

% Windowing
WS=exp(-(0.5*sw^2)*(tt-tt(InS)).^2).*SynC;  % Ventana en la correlacion

% Parameters for "bootstraping"
InSR = InS;

% Compute for all frequencies
for nF=1:length(Freq)
    
    si=0.1*2*pi*Freq(nF);   % Sigma i
    % Crosscorrelagram and Autocorrelagram filtering
    dS=computebandfftfilter_gauss(WS',dt,Freq(nF),si);
    
    % Let's do bootstraping
    [~,InS] = max(dS);
    
    BS = 1;
    while BS == 1
        if (tt(InS) < tt(InSR)+0.5/Freq(nF)) && (tt(InS) >= tt(InSR)-0.5/Freq(nF))
            BS = 0;
        elseif tt(InS) >= tt(InSR)+0.45/Freq(nF)
                InS=InS-round(1/Frec(nF)/dt);
        elseif tt(InS) < tt(InSR)-0.45/Freq(nF)
                InS=InS+round(1/Freq(nF)/dt);
        end
    end
    
    % Model for fitting
    Eqn = @(b,x)(b(1)*exp(-0.5*b(2).^2*(x-b(3)).^2).*cos(b(4)*(x-b(5))));
    [As,Is] = max(envelope(dS));
    GaS = nlinfit(tt,dS,Eqn,[As 2*pi*si tt(Is) 2*pi*Freq(nF) tt(InS)]);
    

    wf(nF) = GaS(4);
    sf(nF) = GaS(2);
end