%%          The GLOBEC Kriging Software Package - EasyKrig

%% Copyright (c) 1998, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%% Institution.  All Rights Reserved.
%% The kriging software described in this document was developed by Dezhang Chu
%% with funding from the National Science Foundation through the U.S. GLOBEC
%% Georges Bank Project's Program Service and Data Management Office.  It was
%% inspired by a MATLAB toolbox in MATLAB developed by Yves Gratton and
%% Caroline Lafleur (INRS-Océanologie, Rimouski, Qc, Canada) and Jeff Runge
%% (Institut Maurice-Lamontagne).  The original version of trans.m was written by D.
%% Marcotte and is used with permission.  Also, we are using, with permission, the
%% variogram model code from 'variogr2.m' written by Yves Gratton and Caroline
%% Lafleur.  This software may be reproduced for noncommercial purposes only.

%% This program is distributed in the hope that it will be useful, but WITHOUT ANY
%% WARRANTY.

%% Contact Dr. Chu at dchu@whoi.edu with enhancements or suggestions for changes.



clear global
close all

global para	 		% all setting and processing parameters
global data  		% input and output data
global hdl			% handles for all windows
global color		% color RGB

%% Matlab Version
Ver=ver;
Version=[];
for i = 1:length(Ver)
    p=strcmp(Ver(i).Name,'MATLAB');
    q=strcmp(Ver(i).Name,'MATLAB Toolbox');
    if p == 1 | q == 1
        Version = floor(str2num(Ver(i).Version(1:3)));
        break
    end
end
if isempty(Version)
    disp(['Can''t determine what version of your Matlab software'])
    break
end

para.Matlab_Version=Version;
initialization6x
hdl.navigator.h0=main_menu3d;
disp_images;

