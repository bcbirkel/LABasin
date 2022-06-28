%% rotate_to_gcp.m
% function to rotate east and north component seismograms to radial and
% transverse component seismograms
%
%% usage:
% [radial,transverse] =
% rotate_to_gcp(east_component_vector,north_component_vector,azimuth)
% This instance returns the radial and transverse components assuming an
% east azimuth of 90 degress and north azimuth of 0 degrees
%
%% optional usage:
% [radial,transverse] =
% rotate_to_gcp(east_component_vector,north_component_vector,azimuth, east_azimuth, north_azimuth)
% This usage allows for misaligned or previously rotated sensors.
%
%% second optional usage:
% [radial,transverse,radial_azimuth,transverse_azimuth] =
% rotate_to_gcp(east_component_vector,north_component_vector,azimuth)
% This usage returns the rotation angles computed
%
%% third optional usage:
% [radial,transverse,radial_azimuth,transverse_azimuth] =
% rotate_to_gcp(east_component_vector,north_component_vector,azimuth,east_azimuth, north_azimuth)
% This is the full available implementation
%
%% c version source code stripped of initializations etc...
%         /* setup the rotation angles */
%        //angle between east component and event
%        delta_east_azimuth = east_azimuth - azimuth;
%
%        //angle between north component and event
%        delta_north_azimuth = north_azimuth - azimuth;
%
%        //translate to radians with short names for easy use
%        der = delta_east_azimuth * PI / 180.0;
%        dnr = delta_north_azimuth * PI / 180.0;
%
%        /* now do the rotation */
%        for (i=0; i<sac_header_e.npts; i++) {
%                signal_t[i] = signal_e[i] * cos(dnr) + signal_n[i] * sin(dnr);
%                signal_r[i] = -1.0 * signal_e[i] * sin(dnr) + signal_n[i] * cos(dnr);
%
%        }
%
%%
function [radial,transverse, delta_east, delta_north] = rotate_to_gcp(east,north,azimuth,varargin)
    % check if the user supplies sensor azimuths
    if nargin == 3
            east_azimuth = 90;
            north_azimuth = 0;
    elseif nargin == 4
            disp('Error, must give both components or assumed 90 and 0!');
    elseif nargin == 5
            east_azimuth = varargin{1}(1);
            north_azimuth = varargin{2}(1);
    else
        disp('Error, improper number of inputs given. see the help file');
    end
    
    % check vector lengths
    east_length = length(east);
    north_length = length(north);
    npts = min(east_length, north_length);
   
    % initialize output vectors
    radial = zeros(1,npts);
    transverse = zeros(1,npts);
    
    % setup the values
    der = (east_azimuth - azimuth) * pi / 180;
    dnr = (north_azimuth - azimuth) * pi / 180;
    
    % only do rotate for the points overlapping
    for i=1:npts
        radial(i) = -1.0 * east(i) * sin(dnr) + north(i) * cos(dnr);
        transverse(i) = east(i) * cos(dnr) + north(i) * sin(dnr);
    end
    
    delta_east = der * 180 / pi;
    delta_north = dnr * 180 / pi;
    
    return
 
