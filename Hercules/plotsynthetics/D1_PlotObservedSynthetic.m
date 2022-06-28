% Before running this code make sure your simulation in Hercules has
% finished. Here you'll compare data and synthetics.
% D1_ObvservedSynthetics.m plots the seismograms from the forward simulation
% in Hercules. This program reads the observed database, reads the
% outputplane with the synthetic seismograms and plots both.
% Written by Alan Juarez (UNAM-USC) Nov, 2016.
% Edited Nov, 2019
% Modified by Brianna Birkel, Mar 2020

clear all; close all; clc

% Loading Earthquakes and Stations for the Simulation
load('/home/scec-00/birkel/hercules_preprocessing/MatlabCodes/DataPreparation/EarthquakeDatabase.mat');


% Specify the desired minimum and maximum % frequencies for waveform filtering. 
% (it is a super good idea these match with the hercules simulation)
FMax = 1/5;
FMin = 1/40;
Mo= 4.0972e+30; % Seismic Moment from Hercules
Dt = 0.25;
T = 0:Dt:200;
ns = 1;           % Number of Source to plot (SourceN)
 
% Loading information of the box simulation
FileID = fopen('/home/scec-00/birkel/hercules_preprocessing/MatlabCodes/AuxFiles/XYtoLonLat.txt');

Zone = fscanf(FileID, '%f', 1); Az = fscanf(FileID, '%f', 1);
X0_lat = fscanf(FileID, '%f', 1); Y0_lon = fscanf(FileID, '%f', 1);
XLim(1) = fscanf(FileID, '%f', 1); XLim(2) = fscanf(FileID, '%f', 1);
YLim(1) = fscanf(FileID, '%f', 1); YLim(2) = fscanf(FileID, '%f', 1);
ZLim(1) = fscanf(FileID, '%f', 1); ZLim(2) = fscanf(FileID, '%f', 1);
Corner(1) = fscanf(FileID, '%f', 1); Corner(2) = fscanf(FileID, '%f', 1);
fclose all; An=abs(Az);

nP=1;
for i=1:2
    for j=1:2
        ROT=[cosd(Az) sind(Az); -sind(Az) cosd(Az)]*[YLim(i); XLim(j)]+[Y0_lon; X0_lat];
        [BoxSim(nP,1),BoxSim(nP,2)] = utm2ll(ROT(1),ROT(2),Zone);
        nP=nP+1;
    end
end
BoxSim(nP,:)=BoxSim(1,:);
BoxSim([3 4],:)=BoxSim([4 3],:);


% Ploting recording at each station. Open plane
Plane=fopen('/auto/scec-00/birkel/hercules_preprocessing/Simulations/Source1/output/planes/planedisplacements.0');

% Header and auxiliary data
Aux=  fread(Plane, 8,'double');
BoxCorners=reshape(Aux,4,2);
LengthY=  fread(Plane, 1,'double');
LengthX=  fread(Plane, 1,'double');
PointsDip = fread(Plane, 1,'int');
PointsStrike    = fread(Plane, 1,'int');
DeltaT  = fread(Plane, 1,'double');
Steps      = fread(Plane, 1,'int');
TimeSteps = round(T(end)/(Steps*DeltaT)+1);

% Reading Components for all stations
SynX = zeros(PointsStrike,TimeSteps);
SynY = zeros(PointsStrike,TimeSteps);
SynZ = zeros(PointsStrike,TimeSteps);

for nM = 1:TimeSteps
    for iE = 1:PointsStrike
        SynX(iE,nM) = fread(Plane, 1,'double');
        SynY(iE,nM) = fread(Plane, 1,'double');
        SynZ(iE,nM) = fread(Plane, 1,'double');
    end
end
fclose all;

%% Rotation componets form LonLat to XY and plot
NF = 1; % Number of figure
NS = 1; % Subplot number

% Correct moment
Factor = ddot(Eqk(ns).M,Eqk(ns).M)^0.5/Mo;
    
SynX = Factor*(SynX);
SynY = Factor*(SynY);
SynZ = Factor*(SynZ);
 
% Plot Seismograms
for nr = 1:5%:10%PointsStrike
    
    % Assign components
    SX = rtrend(SynX(nr,:));
    SY = rtrend(SynY(nr,:));
    SZ = rtrend(SynZ(nr,:));
    
    DN = Eqk(ns).Stn(nr).Channels(1).Data;
    DE = Eqk(ns).Stn(nr).Channels(2).Data;
    DZ = Eqk(ns).Stn(nr).Channels(3).Data;
    
    %Take into account if any component has no record
    if norm(DZ)==0
        SZ = 0*SZ;
    end
    if norm(DN)==0 && norm(DE)==0
        SX = 0*SX;
        SY = 0*SY;
    end
    if norm(DN)==0 || norm(DE)==0
        DXS=SX;
        DYS=SY;
        if norm(DN)==0
            DESR = cosd(An)*DYS-sind(An)*DXS;
            DNSR = 0*(sind(An)*DYS+cosd(An)*DXS);
        elseif norm(DE)==0
            DESR = 0*(cosd(An)*DYS-sind(An)*DXS);
            DNSR = sind(An)*DYS+cosd(An)*DXS;
        end
        SY = cosd(An)*DESR+sind(An)*DNSR;
        SX = -sind(An)*DESR+cosd(An)*DNSR;
        clear DXS DYS DNSR DESR;
    end
    
    %Rotate N and E to X and Y
    DY = cosd(An)*DE+sind(An)*DN;
    DX = -sind(An)*DE+cosd(An)*DN;
    clear DN DE;
    
    %Rotate X,Y to R,T
    AnR = atan2(Eqk(ns).Stn(nr).X-Eqk(ns).X,Eqk(ns).Stn(nr).Y-Eqk(ns).Y)*180/pi;
    Obs(1).Signal = cosd(AnR)*DY+sind(AnR)*DX;
    Obs(2).Signal = -sind(AnR)*DY+cosd(AnR)*DX;
    Obs(3).Signal = DZ;
    
    Syn(1).Signal = cosd(AnR)*SY+sind(AnR)*SX;
    Syn(2).Signal = -sind(AnR)*SY+cosd(AnR)*SX;
    Syn(3).Signal = SZ;
    
    clear DX DY DZ SX SY SZ;
    
    t25 = tukeywin(length(Syn(1).Signal),0.05);
    CMP = ['-R';'-T';'-Z'];
    [Am,In] = max(max(abs([Obs(1).Signal; Obs(2).Signal; Obs(3).Signal])));
    
    figure(NF);
    
    for iP=1:3
        Sint = bp_bu_co(t25'.*Syn(iP).Signal,FMin,FMax,1/Dt,4,2);
        Obse = bp_bu_co(t25'.*Obs(iP).Signal,FMin,FMax,1/Dt,4,2);
        
        subplot(3,3,NS+3*(iP-1))
        plot(T,Obse,'k','linewidth',2), hold on
        plot(T,Sint,'r','linewidth',1.5),
        xlabel 'Time [s]', %legend('Observado','Sintetico'), 
        xlim([0 T(end)]), ylim(Am*[-1 1])
        ylabel('Vel [m/s]')
        text(2.5,0.25*Am,[Eqk(ns).Stn(nr).Name  CMP(iP,:)])
    end
    saveas(gcf, 'SynObs_test.png')
    NS=NS+1;
    if NS>3
        NS=1; NF=NF+1;
%         print(['Seismograms' num2str(NF-1) '.pdf'],'-dpdf','-bestfit')
    end
    clear Syn Obs;
end


