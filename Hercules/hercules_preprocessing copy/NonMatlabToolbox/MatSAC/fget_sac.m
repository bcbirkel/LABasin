function [t,data,SAChdr] = fget_sac(filename,earthquaketime,computer)
%[t,data,SAChdr] = fget_sac(filename)

% read sac into matlab 
% written by Zhigang Peng
% program called
% [head1, head2, head3, data]=sac(filename);
% [SAChdr]=sachdr(head1, head2, head3);

% Updated Mon Jul 30 11:21:24 PDT 2001

if nargin <1, error('ERROR!! No input file name'); end

[head1, head2, head3, data.signal]=sac(filename,computer);
[SAChdr]=sachdr(head1, head2, head3);
% Check what quantity it is

data.type=3; % instrumental response included :)
%if( SAChdr.descrip.idep == 5 ) % acc in nm/s2 then the data has to be scaled
%    data.type=2;
%    data.signal=data.signal*1e-9
%end

%if( (SAChdr.descrip.idep==4) ) % vel in volts
%    data.type=4;
%end

%if( SAChdr.descrip.idep == 3 ) % vel in nm/s
%    data.type=1;
%    data.signal=data.signal*1e-9
%end


% compute the time with zero in the event time
%t = [SAChdr.times.b:SAChdr.times.delta:(SAChdr.data.trcLen-1)*SAChdr.times.delta+SAChdr.times.b]';

% Modified by Leonardo Ramirez-Guzman (usgs golden)
originSignal= SAChdr.event.nzhour*60*60+SAChdr.event.nzmin*60+SAChdr.event.nzsec;
originSignal-earthquaketime;
t = 0:SAChdr.times.delta:(SAChdr.data.trcLen-1)*SAChdr.times.delta;
t=t+originSignal-earthquaketime;

