% Full Waveform Inversion - Code No.D0.
% D0_PlotSynthetictest.m plots the synthetics from the forward simulation 
% in Hercules. This program reads outputplane with the synthetics . 

% Written by Alan Juarez (UNAM-USC) Nov, 2016.
clear all; close all; clc

%Path to the Source directory. In this directory are
% saved the outputs from the Hercules Simulation
PathSrc = '/project/scec_608/birkel/runSim/Simulation/AdjointMov/Phase1/';

%Specify the desired minimum and maximum
% frequencies for waveform filtering. (it is a super good idea these match with
% the hercules simulation)
FMax=1;
FMin=1/20;

%Ploting recording at each station. Open plane
PathPlane=[PathSrc '/output/planes/planedisplacements.0'];
Plane=fopen(PathPlane);

% Header and auxiliary data
Aux=  fread(Plane, 8,'double');
BoxCorners=reshape(Aux,4,2);
LengthY=  fread(Plane, 1,'double');
LengthX=  fread(Plane, 1,'double');
PointsDip = fread(Plane, 1,'int');
PointsStrike    = fread(Plane, 1,'int');
DeltaT  = fread(Plane, 1,'double');
Steps      = fread(Plane, 1,'int');
TimeSteps = round(300/(Steps*DeltaT)+1);

% Reading Components for all stations
for nM=1:TimeSteps
    for iE=1:PointsStrike
        SynS(iE).X(nM) = fread(Plane, 1,'double');
        SynS(iE).Y(nM) = fread(Plane, 1,'double');
        SynS(iE).Z(nM) = fread(Plane, 1,'double');
    end
end
fclose(Plane);

Dt=(Steps*DeltaT);
T=0:Dt:300;

% Rotation componets form LonLat to XY and plot

for nE=1:PointsStrike
    Syn(1).Signal = SynS(nE).X;
    Syn(2).Signal = SynS(nE).Y;
    Syn(3).Signal = SynS(nE).Z; %Remember, Hercules saves vertical downward
    
    figure;
    for iP=1:3
        Sint = Syn(iP).Signal;
        Sint=bp_bu_co(Sint,FMin,FMax,1/Dt,4,1);
        subplot(3,1,iP)
        plot(T,Sint,'linewidth',2),
        xlabel 'Time [s]', grid on
    end
    clear Syn;
end
