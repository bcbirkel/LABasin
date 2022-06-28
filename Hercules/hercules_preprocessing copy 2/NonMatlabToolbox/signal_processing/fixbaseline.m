function [signalfinal]=fixbaseline(signal,p)
% obtain the best p fit
[p,ErrorEst] = polyfit(1:length(signal),signal',p);
x=1:length(signal);
signalfinal=(signal'-polyval(p,x))';

return