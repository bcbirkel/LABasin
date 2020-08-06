function [S_BASE,C_TF] = atom_decm_gauss(ds,dt,FMin,FMax,K)
% [BASE,TF] = decomposition(ds,dt,FMin,FMax,K) estimates a seismogram 
% decomposition into seismic phases that are localized in time and
% frequency.
%
% Inputs:
% ds is the sinthetic seismogram.
% dt is the sampling interval.
% Fmin and Fmax are the allowed frequencies for the decomposition (useful
%               to avoid many tiny phases at short frequencies.
% K is the scaling factor for the S-transform (See Masinha, 1998)
%
% Output:
% BASE is a matrix, each row of BASE contains the m-th isolated phase
%       of the sinthetic seismogram. ds~sum(BASE)
% TF is a m*2 matrix with the centre time and centre frequency index
%       voice in BASE.


ns = length(ds);        % Number of samples
ds = reshape(ds,1,ns);  % Sort data (for sts)
T = 0:dt:dt*(ns-1);     % Time vector for plots
F = K*(0:floor(ns/2))/2/floor(ns/2)/dt;  % Frequency vector

sts = stran(ds,K);      % S-transform of the original signal
Ws = abs(sts).^2;       % Power spectrum

Signif = spec_signif(ds,0.01);  % Spectral significancy at 99% (works good for synthetic data)
Sig95x = (Signif')*(ones(1,ns));% expand signif --> (J+1)x(N) array

fx = 0*Sig95x;          % Significance filter
fx(Ws > Sig95x) = 1;    % It's sharp, maybe smooth

[~,iMax,~,~] = extrema2(fx.*Ws);% Find local maxima in the power spectrum
[locf,loct] = ind2sub(size(fx),iMax);   
[~,In]=sort(loct);      % Sort local maxima according to arrival times
loct=loct(In);
locf=locf(In);

%   Set acceptable frequencuencies (Cone of influence)
for i = 1:length(loct)  
    if T(loct(i)) < 0.5/F(locf(i)) || T(loct(i)) > T(end)-0.5/F(locf(i)) % Arribal times should be later than a period
        loct(i) = 0;
    end                 
    if (F(locf(i)) >= FMax) || (F(locf(i)) <= FMin)     % Freqs have to be in the prefered interval
        loct(i) = 0;
    end    
end
locf(loct==0)=[];       % Delete the points that we don't need
loct(loct==0)=[];
loct(locf>length(F)/2)=[];
locf(locf>length(F)/2)=[];

% Gaussian filter parameters
[X1,X2] = meshgrid(T,F);

r = 1;      % Phase counter (not all the values will be keept)
% figure;
for i = 1:length(loct)
    TF = [T(loct(i)) F(locf(i))];
    Sigma = [1/F(locf(i))/4 0; 0 4*F(locf(i))];
    iFilt = mvnpdf([X1(:) X2(:)],TF,Sigma);
    iFilt = reshape(iFilt,length(F),length(T));
    
%     subplot(2,1,1)
%     surf(T,F,log2(Ws))
%     shading interp, view (2), set(gca,'yscale','log'), hold on
%     ylim([2*F(2) F(end)/2]), xlim([T(1) T(end)])
%     xlabel('Time (s)'), ylabel('Frequency (Hz)'), 
%     CM = loadcmap('BlWhYeRe.c3g'); colormap(CM),
% 
%     subplot(2,1,2)
%     surf(T,F,iFilt)
%     shading interp, view (2), set(gca,'yscale','log'), hold on
%     ylim([2*F(2) F(end)/2]), xlim([T(1) T(end)])
%     xlabel('Time (s)'), ylabel('Frequency (Hz)'), 
%     CM = loadcmap('BlWhYeRe.c3g'); colormap(CM),

    FiltroI = istransform(iFilt.*sts)'; % Estimate inverse transform of isolated phase
    
    if norm(FiltroI) > 0.001*norm(ds)	% If its amplitude is too small just reject
        BASE(r,:) = FiltroI;            % Keep phase in basis
        CTF(r,:) = [T(loct(i)) F(locf(i))]; % Save central time and frequency (just in case)
        r = r+1;        % Increment basis counter
    end
    
end
    
[~,InM] = sort(CTF(:,1));   % Sort arrival times

%Sort according travel times
S_BASE = BASE(InM,:);
C_TF = CTF(InM,:);
% 
figure;
subplot(2,2,1)
plot(T,ds,'linewidth',2)
xlabel 'Time (s)', ylabel 'Amplitude'

subplot(2,2,3)
surf(T,F,log2(Ws))
shading interp, view (2), set(gca,'yscale','log'), hold on
ylim([2*F(2) F(end)/2]), xlim([T(1) T(end)])
xlabel('Time (s)'), ylabel('Frequency (Hz)'), 
CM = loadcmap('BlWhYeRe.c3g'); colormap(CM),

subplot(2,2,3)
stem3(T(loct),F(locf),1e10*ones([1 length(loct)]),'p','markersize',12,...
    'MarkerEdgeColor','k','MarkerFaceColor',0.75*[1 1 1],'linewidth',1)

subplot(2,2,[2 4])
for np = 1:size(S_BASE,1)
    plot(T,np+S_BASE(np,:)/max(S_BASE(:)),'color',[0.8500, 0.3250, 0.0980],'linewidth',2), hold on
end
xlabel 'Time (s)', ylabel 'Phases'

end
