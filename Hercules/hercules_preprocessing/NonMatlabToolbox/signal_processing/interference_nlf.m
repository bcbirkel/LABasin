function [GSDF,Jpq,C,B,P] = interference_nlf(Base,Syn,Obs,dt,Freq)
%
% FWC = interference(Base,Dt,Frec) estimates the interference factors of
% the synthetic waveform dictionary contained in Base. See Gee & Jordan
% 1992 for a detail explanation of the interference.
%
% Input:
%   Base - is a m*n array, where m is the number of elements of the
%       dictionary and n is the number of time samples.
%   Dt - is the time sampling of the time series.
%   Frec - is a array with the frequencies at which the intereference is
%   estimated.
% Output:
%   C - are the interference factors
%   GSDF(m,f,t) - measured GSDF for the m-phase, f-frequency, and t-GSDF
% The total interference can be estimated as C = B.*cos(P);
%
% This function uses the nlinfit function. Results can be innaccurate if
% warnings are shown.

% Data, arrays, and parameters needed
% Syn = sum(Base,1)';
T = dt*(0:length(Syn)-1)';
sw = mean(Freq);	% Effective bandwith

% Estimate GSDF for output
GSDF = zeros([size(Base,1) length(Freq) 5]);
Jpq = zeros([size(Base,1) length(Freq) length(Syn)]);
for nb = 1:size(Base,1)
    Iso = Base(nb,:)';
    Amp = 10/max(abs(Iso)); % Normalization factor (helps to reduce numerical underflow)
    
    [GSDF(nb,:,:),Jpq(nb,:,:),~,~] = gsdf_simple(Amp*Syn,Amp*Obs,Amp*Iso,dt,Freq);
    Jpq(nb,:,:) = Amp*Jpq(nb,:,:);    % Correct normalization
end

% Estimate GSDF for each combination of Isol Filters at all frequencies
for m = 1:size(Base,1)
    for n = 1:size(Base,1)
        Amp = 1/min(max(abs(Base(n,:))),max(abs(Base(m,:))));    % Normalization factor (helps to reduce numerical underflow)
        
        [~,~,FWCmn(m,n).PS,~] = gsdf_simple(Amp*Base(n,:)',Amp*Base(m,:)',Amp*Base(m,:)',dt,Freq);
        FWCmn(m,n).PS(:,1) = FWCmn(m,n).PS(:,1)/Amp/Amp;    % Correct normalization
    end
end

% Estimate a normalization factor (helps to reduce numerical underflow)
Amp = 10/max(abs(Syn));
        
% Effective frequency and bandwidth
[~,~,PS,~] = gsdf_simple(Amp*Syn,Amp*Obs,Amp*Syn,dt,Freq);
PS(:,1) = PS(:,1)/Amp/Amp;

[wf,sf] = compute_wf_sf_nlf(Amp*Syn,dt,Freq);

% Synthetic parameters
tp = PS(:,5);
tg = PS(:,3);
tq = -log(PS(:,1))./wf(:);
ta = -(PS(:,4)-wf(:))./(sf(:).^2);

% Interference parameters
B = zeros([size(Base,1) size(Base,1) length(Freq)]);
P = B;
for m = 1:size(Base,1)
    for n = 1:size(Base,1)
        for nF = 1:length(Freq)
            FWCmn(m,n).tp(nF) = FWCmn(m,n).PS(nF,5);
            FWCmn(m,n).tg(nF) = FWCmn(m,n).PS(nF,3);
            FWCmn(m,n).tq(nF) = -log(FWCmn(m,n).PS(nF,1))/wf(nF);
            FWCmn(m,n).ta(nF) = -(FWCmn(m,n).PS(nF,4)-wf(nF))/sf(nF)^2;
            
            B(m,n,nF) = exp(-wf(nF)*(FWCmn(m,n).tq(nF)-tq(nF)))...
                *exp(-0.5*sf(nF)^2*(FWCmn(m,n).tg(nF)-tg(nF))^2);
            P(m,n,nF) = (wf(nF)-sf(nF)^2*FWCmn(m,n).ta(nF))*...
                (FWCmn(m,n).tp(nF)-tg(nF))-wf(nF)*(tp(nF)-tg(nF));
        end
    end
end

C = B.*cos(P);


