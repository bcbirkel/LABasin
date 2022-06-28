% Full Waveform Inversion - Code No.D1.
% Before running this code make sure your simulation in Hercules has
% finished. Here you'll compare data and synthetics.
% D1_ObvservedSynthetics.m plots the seismograms from the forward simulation
% in Hercules. This program reads the observed database, reads the
% outputplane with the synthetic seismograms and plots both.
% Special care must be taken in the folder structure where catalogs,
% resources and databases are located or where they are going to be saved.
% This is probably a mess, but it's helpfull to keep everything in order.

% Written by Alan Juarez (UNAM-USC) Nov, 2016.
clear all; close all; clc

%% UserModification: Path to the Source directory. In this directory are
% saved the outputs from the Hercules Simulation
PathSrc = '/project/scec_608/birkel/S20140329/Simulation/Forward/';

%% UserModification: Path to the Resources directory (catalogs, stations...)
PathRes='/project/scec_608/birkel/S20140329/Resources/';

%% UserModification: Path to the AuxFiles directory. Here we save or read
%  things that previous codes saved, for example the database.
PathAux = '/project/scec_608/birkel/S20140329/Codes/AuxFiles/';

% Loading Earthquakes and Stations for the Simulation
load([PathAux 'EarthquakeDB.mat']);


%% UserModification: Specify the desired minimum and maximum
% frequencies for waveform filtering. (it is a super good idea these match with
% the hercules simulation)
FMax=1;
FMin=1/100;
Mo= 4.065432e+30; % Seismic Moment from Hercules
Dt=0.25;
T=0:Dt:200;
iS=1;           % Number of Source to plot (SourceN)
NumEst=[47    39    29    24    11     8    59    48    32    56    69    65 ...
        12    51    20    43    25    60    30    28    17    64    18    22 ...
        62    57    50    46    26    58    42    68    34    23    52    63 ...
        31    36    45    37    66    33    55    40    44    67    38    49 ...
        27    61    35    41    19     5    70    53    54     4    15    71];
 
%% No modification needed in this secction 16
% Loading information of the box simulation
FileID = fopen([PathAux 'XYtoLonLat.txt']);

Zone = fscanf(FileID, '%f', 1); Az = fscanf(FileID, '%f', 1);
X0_lat = fscanf(FileID, '%f', 1); Y0_lon = fscanf(FileID, '%f', 1);
XLim(1) = fscanf(FileID, '%f', 1); XLim(2) = fscanf(FileID, '%f', 1);
YLim(1) = fscanf(FileID, '%f', 1); YLim(2) = fscanf(FileID, '%f', 1);
ZLim(1) = fscanf(FileID, '%f', 1); ZLim(2) = fscanf(FileID, '%f', 1);
Corner(1) = fscanf(FileID, '%f', 1); Corner(2) = fscanf(FileID, '%f', 1);
fclose(FileID); An=abs(Az);

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


%% No modification needed in this secction
% Ploting recording at each station. Open plane
PathPlane=[PathSrc '/Source1/output/planes/planedisplacements.0'];
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
TimeSteps = round(T(end)/(Steps*DeltaT)+1);

% Reading Components for all stations
for nM=1:TimeSteps
    for iE=1:PointsStrike
        SynS(iE).X(nM) = fread(Plane, 1,'double');
        SynS(iE).Y(nM) = fread(Plane, 1,'double');
        SynS(iE).Z(nM) = fread(Plane, 1,'double');
    end
end
fclose(Plane);

% Rotation componets form LonLat to XY and plot
NF=1; % Number of figure
NS=1; %Subplot number
cd Figures
for nE=NumEst
    if Mo == 0
        Factor = 1;
    else
        Factor = (Eqk(iS).Src.Mo*10^Eqk(iS).Src.Exp)/Mo;
    end
    SX = Factor*(SynS(nE).X);
    SY = Factor*(SynS(nE).Y);
    SZ = -Factor*(SynS(nE).Z);

    if length(Eqk(iS).Stn(nE).Cmpt) == 3
        DN=Eqk(iS).Stn(nE).Cmpt(1).Signal;
        DE=Eqk(iS).Stn(nE).Cmpt(2).Signal;
        DZ=Eqk(iS).Stn(nE).Cmpt(3).Signal;
    elseif length(Eqk(iS).Stn(nE).Cmpt) == 2
        DN=Eqk(iS).Stn(nE).Cmpt(1).Signal;
        DE=Eqk(iS).Stn(nE).Cmpt(2).Signal;	
	DZ=0*DN;
    elseif length(Eqk(iS).Stn(nE).Cmpt) == 1
	DN=Eqk(iS).Stn(nE).Cmpt(1).Signal;
	DE=0*DN;
	DZ=0*DN;
    end

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
    DY=cosd(An)*DE+sind(An)*DN;
    DX=-sind(An)*DE+cosd(An)*DN;
    clear DN DE;
    
    %Rotate X,Y to R,T
    AnR=atan2(Eqk(iS).Stn(nE).X-Eqk(iS).Loc.X,Eqk(iS).Stn(nE).Y-Eqk(iS).Loc.Y)*180/pi;
    Obs(1).Signal = cosd(AnR)*DY+sind(AnR)*DX;
    Obs(2).Signal = -sind(AnR)*DY+cosd(AnR)*DX;
    Obs(3).Signal = DZ;
    
    Syn(1).Signal = cosd(AnR)*SY+sind(AnR)*SX;
    Syn(2).Signal = -sind(AnR)*SY+cosd(AnR)*SX;
    Syn(3).Signal = SZ;
    
    clear DX DY DZ SX SY SZ;
    
    t25 = tukeywin(length(Syn(1).Signal),0.05);
    CMP = ['-R';'-T';'-Z'];
    [Am,In]=max(max(abs([Obs(1).Signal; Obs(2).Signal; Obs(3).Signal])));
    t0=max([40 min([5*floor(T(In)/5) 200-80])]);
    
    figure(NF);
    Error(nE)=0;
    Eqk(iS).Stn(nE)
    for iP=1:3
        Sint=bp_bu_co(t25'.*Syn(iP).Signal,FMin,FMax,1/Dt,4,2);
        Obse=bp_bu_co(t25'.*Obs(iP).Signal,FMin,FMax,1/Dt,4,2);
        
        subplot(3,3,NS+3*(iP-1))
        plot(T,Obse,'k','linewidth',2), hold on
        plot(T,Sint,'r','linewidth',1.5),
        xlabel 'Time [s]', %legend('Observado','Sintetico'), 
        xlim([t0-40 t0+80]), ylim(Am*[-1 1])
        ylabel('Vel [m/s]')
        text(t0+2.5-40,0.25*Am,[Eqk(iS).Stn(nE).Name  CMP(iP,:)])
        Error(nE)=Error(nE)+(Sint-Obse)*(Sint-Obse)'/(Sint*Sint') ...
            +(Sint-Obse)*(Sint-Obse)'/(Obse*Obse');
    end
    NS=NS+1;
    if NS>3
        NS=1; NF=NF+1;
         print(['Seismograms' num2str(NF-1) '.pdf'],'-dpdf','-bestfit')
    end
    clear Syn Obs;
end


%% No modification needed in this secction
% loading topography file
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

%% Modifications needed in this section.
% Map settings. Just to plot a map with stations that recorded.
[Lon,Lat]=meshgrid(lon,lat);    % lon and lat come from elevationdataforaxis.mat
% We create a colormap that looks nice for topography but its optional
% Creating and setting Map parameters. Here we are using information in the
% resources directory
figure;
ColorMap=[loadcmap('tpglpom.c3g'); loadcmap('GrayWhite.c3g')];
m_proj('Miller','lon',[min(BoxSim(:,2))-0.5 max(BoxSim(:,2))+0.5],...
    'lat',[min(BoxSim(:,1))-0.5 max(BoxSim(:,1))+0.5]);
hold on;

% Plot topography and add some other nice things
m_surf(Lon,Lat,bathymetry,bathymetry),
view (2), colormap(ColorMap),
caxis([-1.1*max(max(bathymetry)) 1.1*max(max(bathymetry))]), shading interp
m_contour(lon,lat,bathymetry,[0, 0],'Color',[0.2 0.2 0.2],'linewidth',1);

%% User modification needed: Choose and edit paths for political limits or
% Fisiographycprovinces. I'm sure this will give some troubles when using
% different regions. But if you know what to do whit these sample files you
% are done.

load([PathRes 'RollinsData/caldata.mat']);
m_plot(uclafaults(:,1),uclafaults(:,2),'color',[0.5 0.5 0.5],'linewidth',1)
% m_plot(stateborders(:,2),stateborders(:,1),'color',[0.1 0.1 0.1],'linewidth',2)


%% No Modifications needed in this section. Plot simulation box
m_grid('tickdir','in','box','fancy','fontsize',12,'fontweight','demi');
m_plot(BoxSim(:,2),BoxSim(:,1),'color',[0.7 0.1 0.1],'linewidth',3)


%% No Modifications needed in this section.
% Here we are going to pick all sources and stations inside the simulation
% box. It is probably that something wont look fine when plotting map
CS=[.64 .08 .18];
m_plot(Eqk(iS).Loc.Lon,Eqk(iS).Loc.Lat,'p','markerfacecolor',...
    CS,'markeredgecolor',[0 0 0],...
    'MarkerSize',14,'color','k'); hold on

iSt=1; CV=[.47 .67 .19];
for StN=NumEst%1:size(Eqk(iS).Stn,2)
    m_plot(Eqk(iS).Stn(StN).Lon,Eqk(iS).Stn(StN).Lat,'v',...
        'markerfacecolor',CV,'markeredgecolor',[0 0 0],...
        'MarkerSize',8,'color','k'); hold on
    m_text(Eqk(iS).Stn(StN).Lon+0.2,Eqk(iS).Stn(StN).Lat,[Eqk(iS).Stn(StN).Name])
end

% print(['MapStnEqk.pdf'],'-dpdf','-bestfit')
cd ..

% %% No Modifications needed in this section.
% % Plot legends and some other fancy things to the MAP
% MinLat = min(BoxSim(:,1))-0.5;
% MinLon = min(BoxSim(:,2))-0.5;
% m_plot(MinLon+0.2,MinLat+0.5,'p',...
%     'markerfacecolor',CS,'markeredgecolor',[0 0 0],...
%     'MarkerSize',10,'color','k'); hold on
% m_plot(MinLon+0.2,MinLat+0.2,'v',...
%     'markerfacecolor',CV,'markeredgecolor',[0 0 0],...
%     'MarkerSize',8,'color','k'); hold on
% 
% m_text(MinLon+0.2,MinLat+0.5,'  Epicenter','FontSize',10)
% m_text(MinLon+0.2,MinLat+0.2,'  Receivers','FontSize',10)
