%
%
% [itphase, dt]=pickfirstarrival(a,dt)
%
%  a  =  signal
%  dt = delta time
%  ll = lower limt of total energy (0-1)
%  ul = upper limt of total energy (0-1)
function [Dk,Ftkk,pickindex,picktime]=pickfirstarrival(signal,dt,originTime)

numSamples=length(signal);
a=signal.*signal;
ck=cumsum(a);
cT=ck(length(ck));
T=(length(a)-1)*dt;

Dk(1)=0;
Ftkk(1)=0;

k=(0:dt:(numSamples-1)*dt);
F1=(ck/cT); F2=(k/T);
Dk=F1-F2;
Ftkk=((cT-ck)./(T-k))./(ck./k);
ckI=find(ck==0);
Ftkk(ckI)=0;
    

[val,pickindex]=max(abs(Dk.*Ftkk));
picktime=originTime+(pickindex-1)*dt;

% figure ()
%   subplot (3,1,1)
%   plot(k,F1,'linewidth',2), xlim([100 500])
%   subplot (3,1,2)
%   plot(k,Dk,'linewidth',2), xlim([100 500])

end
