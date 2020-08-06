% Full Waveform Inversion - Code No.D0.
% D0_PlotSynthetictest.m plots the synthetics from the forward simulation 
% in Hercules. This program reads outputplane with the synthetics . 

% Written by Alan Juarez (UNAM-USC) Nov, 2016.
% Modified by Brianna Birkel, Mar 2020
clear; close all; clc

%Path to the Source directory. In this directory are
% saved the outputs from the Hercules Simulation
PathSrc1 = '/Users/bcbirkel/Documents/Seismofiles/hercules_preprocessing/Simulations/Source1/';
PathSrc2 = '/Users/bcbirkel/Documents/Seismofiles/hercules_preprocessing/Simulations/Source1/';

%Specify the desired minimum and maximum
% frequencies for waveform filtering. (it is a super good idea these match with
% the hercules simulation)
FMax=1/4;
FMin=1/30;

% Load data Source1
PathPlane=[PathSrc1 '/output/planes/planedisplacements.0'];
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
TimeSteps = round(150/(Steps*DeltaT)+1);

% Reading Components for all stations
for nM=1:TimeSteps
    for iE=1:PointsStrike
        SynS1(iE).X(nM) = fread(Plane, 1,'double');
        SynS1(iE).Y(nM) = fread(Plane, 1,'double');
        SynS1(iE).Z(nM) = fread(Plane, 1,'double');
    end
end
fclose(Plane);


%{
% Load data Source2
PathPlane = [PathSrc2 '/output/planes/planedisplacements.0'];
Plane=fopen(PathPlane);

% Header and auxiliary data
Aux= fread(Plane, 8,'double');
BoxCorners=reshape(Aux,4,2);
LengthY=  fread(Plane, 1,'double');
LengthX=  fread(Plane, 1,'double');
PointsDip = fread(Plane, 1,'int');
PointsStrike    = fread(Plane, 1,'int');
DeltaT  = fread(Plane, 1,'double');
Steps      = fread(Plane, 1,'int');
TimeSteps = round(235/(Steps*DeltaT)+1);

% Reading Components for all stations
for nM=1:TimeSteps
    for iE=1:PointsStrike
        SynS2(iE).X(nM) = fread(Plane, 1,'double');
        SynS2(iE).Y(nM) = fread(Plane, 1,'double');
        SynS2(iE).Z(nM) = fread(Plane, 1,'double');
    end
end
fclose(Plane);
%}

% Some simulation parameters
Dt=(Steps*DeltaT);
T=0:Dt:235;
Mo= 2.804148e+30; % Seismic Moment from Hercules
load('/Users/bcbirkel/Documents/Seismofiles/hercules_preprocessing/MatlabCodes/AuxFiles/EqkStnDB.mat');
%Factor = (Eqk(1).Src.Mo*10^Eqk(1).Src.Exp)/Mo;
Factor = 1 %this is just some crap because i can't find EqkStnDB
t25 = tukeywin(length(T),0.25);
% SR = differentiate(triangularPulse(0,((Eqk.Src.Mo*10^(Eqk.Src.Exp))/(2.2e16))^0.333,T));
% SR = fft(SR/max(SR));          % receiver source function
% SR = (conj(SR)./(1e-0*max(abs(SR))^2+SR.*conj(SR)));

iS = 1:size(Eqk.Stn,2)

iS = 1:63
%iS = 1:size(SynS1(1).Z)

% Rotation componets form LonLat to XY and plot
for r = iS
    
%     SynS1(r).Z = ifft(SR.*fft(SynS1(r).Z));
%     SynS2(r).Z = ifft(SR.*fft(SynS2(r).Z));

    Signal1 = Factor*bp_bu_co(SynS1(r).Z,FMin,FMax,1/Dt,4,2);
%    Signal2 = Factor*bp_bu_co(SynS2(r).Z,FMin,FMax,1/Dt,4,2);
    
    Dist(r) = 0.001*norm([Eqk.X Eqk.Y]-[Eqk.Stn(r).X Eqk.Stn(r).Y])
    %Dist_row=size(Dist(r),1)
    %Dist_col=size(Dist(r),2)
    %T_row=size(T,1)
    %T_col=size(T,2)
    %T(r)
    %Dist(r)
    plot(T(r),Dist(r)+4*Signal1/max(Signal1),'color',[0 0 0],'linewidth',1.5), hold on,
%    plot(T,Dist(r)+4*Signal2/max(Signal1),'color',[0.635 0.078 0.184], 'linewidth',1),    
end
xlabel 'Time (s)', ylabel 'Distance (km)', grid on
ylim([0 370]), xlim([0 200])


%% UserModification: Path to the Resources directory (catalogs, stations...)
PathRes='/Users/bcbirkel/Documents/Seismofiles/hercules_dependencies/Resources/';

%% UserModification: Path to the AuxFiles directory. Here we save or read
%  things that previous codes saved.
PAuxFiles = '/Users/bcbirkel/Documents/Seismofiles/hercules_preprocessing/MatlabCodes/AuxFiles/';


%% No modification needed in this secction
% Loading information of the box simulation
FileID = fopen([PAuxFiles 'XYtoLonLat.txt']);

Zone = fscanf(FileID, '%f', 1); Az = fscanf(FileID, '%f', 1);
X0_lat = fscanf(FileID, '%f', 1); Y0_lon = fscanf(FileID, '%f', 1);
XLim(1) = fscanf(FileID, '%f', 1); XLim(2) = fscanf(FileID, '%f', 1);
YLim(1) = fscanf(FileID, '%f', 1); YLim(2) = fscanf(FileID, '%f', 1);
ZLim(1) = fscanf(FileID, '%f', 1); ZLim(2) = fscanf(FileID, '%f', 1);
Corner(1) = fscanf(FileID, '%f', 1); Corner(2) = fscanf(FileID, '%f', 1);
fclose(FileID);

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

%{
%% UserModification: Path to load topography
% In my case I have fancy files with topography, you can edit this
% seccion according to your needs. Topography can be easyly included by 
% downloading topo files from ETOPO1.

% Topography
% loading topography file (if got it from ETOPO you are done!)
elev=load ([PathRes 'Topography/etopo1_bedrock.xyz']);
lon = sort(unique(elev(:,1)));
lat = sort(unique(elev(:,2)),'descend');
i=1;
for jj=1:length(lat)
    for ii=1:length(lon)
        bathymetry(jj,ii)=elev(i,3);
        i=i+1;
    end
end


%% No Modifications needed in this section.
% Map setting
[Lon,Lat]=meshgrid(lon,lat);    % lon and lat come from elevationdataforaxis.mat
% We create a colormap that looks nice for topography but its optional
% Creating and setting Map parameters. Here we are using information in the
% resources directory
figure;
m_proj('Miller','long',[-117 -114],'lat',[32 35]); hold on;
m_surf(Lon,Lat,bathymetry+30,bathymetry+30), shading interp
ColorMap=createcmap([140 200 245]/250,[0.5 0.5 0.5]);
colormap(ColorMap), caxis([-1 1])
m_contour(lon,lat,bathymetry+30,[0, 0],'Color',[0.2 0.2 0.2],'linewidth',1);
m_grid('box','fancy','tickdir','in','fontsize',9,'fontweight','normal'); 
% loading topography file
PathRes='/home/scec-00/alanjuar/hercules_dependencies/Resources/';

elev=load ([PathRes 'Topography/etopo1_bedrock.xyz']);
lon = sort(unique(elev(:,1)));
lat = sort(unique(elev(:,2)),'descend');
i=1;
for jj=1:length(lat)
    for ii=1:length(lon)
        bathymetry(jj,ii)=elev(i,3);
        i=i+1;
    end
end

m_contour(lon,lat,bathymetry+30,[0, 0],'Color',[0.3 0.3 0.3],'linewidth',1);

load([PathRes 'RollinsData/caldata.mat']);
m_plot(uclafaults(:,1),uclafaults(:,2),'color',[0.5 0.5 0.5],'linewidth',2)

bd = load('borderdata.mat');
lat = bd.lat([247:273 276 277 279:282 284:285 287:299 302]); 
lon = bd.lon([247:273 276 277 279:282 284:285 287:299 302]); 
for i =1:size(lat,1)
    m_plot(lon{i},lat{i},'color',[0.1 0.1 0.1],'linewidth',2)
end
%}

% Plot Eqk and stations
for r = iS
    m_plot(Eqk.Stn(r).Lon,Eqk.Stn(r).Lat,'v',...
        'markerfacecolor',[0 .45 .74],'markeredgecolor',[0 0 0],...
        'MarkerSize',7,'color','k'); hold on
%     m_text(Stn(r).Lon+0.05,Stn(r).Lat,num2str(r))
end

for r = [34 37 111 11 100 106 74 57]
    m_plot(Eqk.Stn(r).Lon,Eqk.Stn(r).Lat,'v',...
        'markerfacecolor',[0.4660, 0.6740, 0.1880],'markeredgecolor',[0 0 0],...
        'MarkerSize',7,'color','k'); hold on
    m_text(Eqk.Stn(r).Lon+0.05,Eqk.Stn(r).Lat,[Eqk.Stn(r).Name])
end

saveas(gcf,'Mid_distanceplot.png')

for s=1:size(Eqk,2)
    m_plot(Eqk(s).Lon,Eqk(s).Lat,'p','markerfacecolor',...
        [.64 .08 .18],'markeredgecolor',[0 0 0],...
        'MarkerSize',12,'color','k'); hold on
end

saveas(gcf,'DistancePlot.png')

