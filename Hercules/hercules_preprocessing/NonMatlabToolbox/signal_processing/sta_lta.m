function t0=sta_lta(signal,dt,BE,STA,LTA,TR)
%
% Triggering algorithm based on the ratio of the Short-Term
% Average absolute value to the Long-Term Average absolute value
% of a detrended signal. Everything is in standard units.
%
% INPUT:
%
% signal   Vector containing the signal
% dt       Sampling interval (s)
% BE       Beginning and end time of signal ([s s])
% STA      Short-term averaging window length (s)
% LTA      Long-term averaging window length (s)
% TR       Value of STA/LTA ratio that triggers
% PNL      Minimum window length of any triggered section (s)

% OUTPUT:
%
% trigt   Matrix with begin and end times of triggered sections (s)
% stav    Short-term average of absolute values of detrended signal
% ltav    Long-term average of absolute values of detrended signal
% ratio   Ratio of short-term to long-term average
% tim1    Time axis for 'ratio'
% tim2    Time axis for 'stav'
% tim3    Time axis for 'ltav'
% trigs   Vector with triggering points, in samples
%
  
 
%figure out how many samples the windows encompass
STAsmp=ceil(STA/dt);
LTAsmp=ceil(LTA/dt);
NPTS  =length(signal);

% Detrend the signal so DC value and trend don't play
signal=detrend(signal);
signal=signal/max(abs(signal));


% Calculate long-term and short-term average and their ratio
  ltaVector=ones(LTAsmp,1);
  staVector=ltaVector*0;
  staVector(1:STAsmp)=1;

  ltaVector=ltaVector./(sum(ltaVector));
  staVector=staVector./(sum(staVector));
  % moving average
  stavC=conv(abs(signal),staVector)*100;
  ltavC=conv(abs(signal),ltaVector)*100;
  
  stav=round(stavC(1:length(signal)));
  ltav=round(ltavC(1:length(signal)));
  
  ltav(1:length(LTAsmp))=(ltavC(LTAsmp));
  stav(1:length(STAsmp))=(stavC(STAsmp));
 
  ratio=stav./ltav;
  iN=isnan(ratio);
  ratio(iN)=0;
  clear iN
  iN=isinf(ratio);
  ratio(iN)=0;
  ratio(1:STA*round(1/dt))=0;
  
  t0=-1;
  while t0<STA
      k=find(ratio>TR,1);
      if isempty(k)==1
          TR=TR-0.1;
          t0=-1;
      else
          t0=k*dt+BE-STA ;
          if t0<=STA
              ratio(k)=0;
              t0=-1;
          end
      end
  end
%   figure ()
%   subplot (3,1,1)
%   plot((0:NPTS-1)*dt,ratio,'linewidth',2), xlim([100 500])
  
end