function [p, P, NHNM, NLNM] = NoiseCurves(flag)

% Low and High Noise Curves
% Reference: USGS Series Open-File Report, 93-322, "Title Observations 
% and modeling of seismic background noise, Author: Peterson, J. R., 1993
% Webpage: http://pubs.er.usgs.gov/usgspubs/ofr/ofr93322
% 
% flag: 0 for acceleration curves and 1 for velocity curves
%
% Ericka Alinne Hernandez-Solano y Xyoli Perez-Campos
% November 8th, 2006


% Low Noise Model
lnm = [ 0.10     -162.36     5.64;
        0.17     -166.7      0.0;
        0.40     -170        -8.3;
        0.80     -166.4      28.9;
        1.24     -168.6      52.48;
        2.40     -159.98     29.81;
        4.30     -141.1      0.0;
        5.0      -71.36      -99.77;
        6.0      -97.26      -66.49;
       10.0     -132.18     -31.57;
       12.0     -205.27     36.16;
       15.60    -37.65      -104.33;
        21.90   -114.37     -47.1;
        31.60   -160.58     -16.28;
        45.0    -187.5      0.0;
        70.0    -216.47     15.7;
        101.0   -185.0      0.0;
        154     -168.34     -7.61;
        328     -217.43     11.9;
        600     -258.28     26.6;
        10000   -346.88     48.75;
        100000   -346.88     48.75];
p = lnm(:,1);
a = lnm(:,2);
b = lnm(:,3);
NLNMacc = a + (b .* log10(p)); 

% High Noise Model
hnm = [ 0.1     -108.73     -17.23;
        0.22    -150.34     -80.50;
        0.32    -122.31     -23.87;
        0.8     -116.85     32.51;
        3.8     -108.48     18.08;
        4.6     -74.66      -32.95;
        6.3     0.66        -127.18;
        7.9     -93.37      -22.42;
        15.4    73.54       -162.98;
        20.0    -151.52     10.01;
        354.8   -206.66     31.63;
        100000  -206.66     31.63];
P = hnm(:,1);
A = hnm(:,2);
B = hnm(:,3);
NHNMacc = A + (B .* log10(P));

% semilogx(p,NLNMacc,P,NHNMacc), xlim ([0.01 179]), ylim([-200 -50])
switch flag
    case 0
        NLNM = NLNMacc;
        NHNM = NHNMacc;
        %hold on, semilogx(p,NLNM,P,NHNM,'linewidth',2), 
    case 1
        NLNM = NLNMacc + 20 * log10(p/(2*pi));  % Velocity
        NHNM = NHNMacc + 20 * log10(P/(2*pi)); 
        %hold on, semilogx(p,NLNM,P,NHNM,'linewidth',2), 
end

