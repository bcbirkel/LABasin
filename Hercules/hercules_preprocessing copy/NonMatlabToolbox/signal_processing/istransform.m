function out = istransform(S)
% USAGE: s=inv_s_tran(in)
%
% function to take the inverse S-transform, assuming that the S-transform
% was obtained using the function stran.m.  This is part of the ST IST
% pair.  For a more complete reference, see IEEE TRANSACTIONS ON SIGNAL
% PROCESSING, VOL. 55, NO. 10, OCTOBER 2007
% 
% 
% NOTE: this code is designed to work with the S-transform code found on
% Matlab central:
% http://www.mathworks.com/matlabcentral/fileexchange/45848-stockwell-transform--s-transform-
%
% =====================================================================
%  INPUTS
% 
%   S   the input S transform, as given by stran.m 
%
% =====================================================================
%  OUTPUTS
%
%  out  a time domain signal 
%
% ......................................................................
%  by Christian Poppeliers
%     East Carolina University
%     Dept. of Geosciences
%  modifications by Alan Juarez
%     University of Souther California

nhaf = fix(size(S,2)/2);

% perform the inner sum
Suma = sum(S(2:end,:),2)';

% fix the sign on the imaginary part
if nhaf*2 == size(S,2);
    ss = [Suma flip(conj(Suma))];
else
    ss = [0 Suma flip(conj(Suma))];
end

% take the ifft to get the time-domain signal
out = real(ifft(ss));



