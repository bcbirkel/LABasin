%WTC(t,y)  Matlab script for WAVELET SPECTRUM of the time series Y(t), 
%
% See "http://paos.colorado.edu/research/wavelets/"
% Written January 1998 by C. Torrence
%
% Modified Oct 1999, changed Global Wavelet Spectrum (GWS) to be sideways,
%   changed all "log" to "log2", changed logarithmic axis on GWS to
%   a normal axis.

function [T,G] = wavelet_ge(t,sst)   % input SST time series

%------------------------------------------------------ Computation

% normalize by standard deviation (not necessary, but makes it easier
% to compare with plot on Interactive Wavelet page, at
% "http://paos.colorado.edu/research/wavelets/plot/"


n = length(sst);
dt = (t(end)-t(1))/(n-1);
time = [0:n-1]'*dt + t(1) ;  % construct time array

variance = std(sst)^2;
sst = (sst - mean(sst))/sqrt(variance) ;

pad = 1;      % pad the time series with zeroes (recommended)
dj = 1/12;    % this will do 12 sub-octaves per octave
s0 = 2*dt;    % this says start at a scale of 2 months
MaxScale=(n*.17)*2*dt; % automaxscale
j1 = round(log2(MaxScale/s0)/dj);    % this says do 7 powers-of-two with dj sub-octaves each
lag1 = 0.72;%ar1(sst);  % 0.72 lag-1 autocorrelation for red noise background
mother = 'Morlet';

% Wavelet transform:
[wave,period,scale,coi] = wavelet(sst,dt,pad,dj,s0,j1,mother);
power = (abs(wave)).^2 ;        % compute wavelet power spectrum

% Significance levels: (variance=1 for the normalized SST)
[signif,fft_theor] = wave_signif(1.0,dt,scale,0,lag1,-1,-1,mother);
sig95 = (signif')*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig95 = power ./ sig95;         % where ratio > 1, power is significant

% Global wavelet spectrum & significance levels:
global_ws = dt*dt*variance*(sum(power')/n);   % time-average over all times
dof = n - scale;  % the -scale corrects for padding at edges
global_signif = wave_signif(variance,dt,scale,1,lag1,-1,dof,mother);

% Scale-average between El Nino periods of 2--8 years
avg = find((scale >= 2) & (scale < 8));
Cdelta = 0.776;   % this is for the MORLET wavelet
scale_avg = (scale')*(ones(1,n));  % expand scale --> (J+1)x(N) array
scale_avg = power ./ scale_avg;   % [Eqn(24)]
scale_avg = variance*dj*dt/Cdelta*sum(scale_avg(avg,:));   % [Eqn(24)]
scaleavg_signif = wave_signif(variance,dt,scale,2,lag1,-1,[2,7.9],mother);

T = period;
G = global_ws;
% % whos
% 
% %------------------------------------------------------ Plotting
% 
% %--- Plot time series
% subplot('position',[0.1 0.75 0.65 0.2])
% plot(time,sst,'r','linewidth',2)
% set(gca,'XLim',xlim(:))
% xlabel('Time (ut)')
% ylabel('S(t) (u)'), grid
% title('a) Time Series')
% hold off
% 
% %--- Plot wavelet power spectrum
% subplot('position',[0.1 0.10 0.65 0.6])
% clims=[0.95*min(min(log2(abs(power/variance)))) 1.05*max(max(log2(abs(power/variance))))];
% imagesc(time,log2(period),log2(abs(power/variance)),clims);  %*** uncomment for 'image' plot
% grid on,
% set(gca,'XLim',xlim(:))
% Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
% set(gca,'YLim',log2([min(period),max(period)]), ...
% 	'YDir','reverse', ...
% 	'YTick',log2(Yticks(:)), ...
% 	'YTickLabel',Yticks, ...
%     'layer','top')
% set(gca,'box','on','layer','top');
% % 95% significance contour, levels at -99 (fake) and 1 (95% signif)
% hold on
% xlabel('Time (ut)')
% ylabel('Period (ut)')
% title('b) Wavelet Power Spectrum')
% [c,h] = contour(time,log2(period),sig95,[-99 1],'k');
% set(h,'linewidth',1.5)
% tt=[time([1 1])-dt*.5;time;time([end end])+dt*.5];
% hcoi=fill(tt,log2([period([end 1]) coi period([1 end])]),'w');
% set(hcoi,'alphadatamapping','direct','facealpha',.5)
% hold off
% 
% %--- Plot global wavelet spectrum
% subplot('position',[0.77 0.10 0.2 0.6])
% plot(global_ws,log2(period),'k')
% hold on
% plot(global_signif,log2(period),'--r'), grid
% xlabel('Power (u^2)')
% title('c) Global Wavelet Spectrum')
% set(gca,'YLim',log2([min(period),max(period)]), ...
% 	'YDir','reverse', ...
% 	'YTick',log2(Yticks(:)), ...
% 	'YTickLabel','')
% set(gca,'XLim',[0,1.25*max(global_ws)])
% hold off
% 
% % for i=1:length(period)
% %     Var2Save{i,:}=[num2str(period(i),'%3.4f'),'    ',num2str(global_ws(i),'%3.4f')];
% % end
% % dlmcell(filename,Var2Save,'delimiter',' ')
           
end



