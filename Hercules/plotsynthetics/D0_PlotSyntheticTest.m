% This program reads outputplane with the synthetics . 

% Written by Alan Juarez (UNAM-USC) Nov, 2018. Edit by Brianna Birkel, Mar 2020.

clear; close all; clc

%Paths to the outputs from the Hercules Simulations
%PathSrc1 = '/home/scec-00/birkel/hercules_preprocessing/Simulations/Source1/';
%PathSrc2 = '/home/scec-00/birkel/rawIRISdata/20140329/SCEDC/';

%local
PathSrc1 = '/Users/bcbirkel/Documents/Seismofiles/hercules_preprocessing/Simulations/Source1/';
PathSrc2 = '/Users/bcbirkel/Documents/Research/LABasin/rawIRISdata/20200422/';

%Specify the desired minimum and maximum frequencies for waveform filtering
FMax=1;
FMin=1/40;

%% Simulation

% Load data from Simulation1
Plane1 = fopen([PathSrc1 '/output/planes/planedisplacements.0']);

Aux=  fread(Plane1, 8,'double')
BoxCorners=reshape(Aux,4,2)
LengthY=  fread(Plane1, 1,'double')
LengthX=  fread(Plane1, 1,'double')
PointsDip = fread(Plane1, 1,'int')
PointsStrike    = fread(Plane1, 1,'int')
DeltaT  = fread(Plane1, 1,'double')
Steps      = fread(Plane1, 1,'int')
TimeSteps = round(200/(Steps*DeltaT)+1)
ncount=0;
rcount=0;

for n = 1:TimeSteps
    ncount = ncount + 1;
    niter = n;
    for r = 1:PointsStrike
        rcount = rcount + 1;
        riter = r;
        SynS1(r).X(n) = fread(Plane1, 1,'double');
        SynS1(r).Y(n) = fread(Plane1, 1,'double');
        SynS1(r).Z(n) = fread(Plane1, 1,'double');
    end
end
fclose(Plane1);

% Simulation parameters
Dt = (Steps*DeltaT);
T = 0:Dt:200;
Mo = 3.553074e+30; % Seismic Moment from Hercules
load(['/Users/bcbirkel/Documents/Seismofiles/hercules_preprocessing/MatlabCodes/DataPreparation/EarthquakeDatabase.mat']);
%Factor = (Eqk(1).Src.Mo*10^Eqk(1).Src.Exp)/Mo;
%Factor = 10^-8;
MoFactor = 1;

%% Observed Data

OD_Factor = 100; %scaling plot to make more readable
samples_per_sec = 100;
Dt = 1/samples_per_sec;

%Load observed data
files=dir([PathSrc2 'CE*'])
for r = size(files, 1)
    filename = files(r).name(1:8)
    PZF = dir([PathSrc2 '/NCEDC/SACPZ*' filename '*']);
    PZFiles(3r:3r+2) = PZF
end

for ii=1:3:size(files,1)
	files(ii)
	files(ii).name
    a = (ii+2)/3;
    
    %info about obs data
    hd = rdSacHead([PathSrc2 files(ii).name]);
    SacHead(a).hd = hd;
    SacHead(a).name = files(ii).name(4:6);
    
    [data,H]=rdSac([PathSrc2 files(ii).name]);
    %data = detrend(data(2,:));
    ObsData(a).E = data.';
    ObsDataHead(a).E = H;
    PZ(a).E = PZFiles(ii).name;

    [data,H]=rdSac([PathSrc2 files(ii+1).name]);
    %data = detrend(data(2,:));
    ObsData(a).N = data.';
    PZ(a).N = PZFiles(ii+1).name;

    [data,H]=rdSac([PathSrc2 files(ii+2).name]);
    %data = detrend(data(2,:));
    ObsData(a).Z = data.';
    PZ(a).Z = PZFiles(ii+2).name;
    
	%plotdata(0:a-1,data);
	%saveas(gcf,['SacObs-' ii '.png'])
end

 %{
for r=1:size(ObsData,2)
    % Correct instrument response
    ObsData(r).N = transfer(ObsData(r).N,length(ObsData(r).N),samples_per_sec,PZ(r).N,FMin,FMax);
    ObsData(r).E = transfer(ObsData(r).E,length(ObsData(r).E),samples_per_sec,PZ(r).E,FMin,FMax);
    ObsData(r).Z = transfer(ObsData(r).Z,length(ObsData(r).Z),samples_per_sec,PZ(r).Z,FMin,FMax);

    %ObsData(r).N = [diff(cObsData(r).N)]/Dt;
    %ObsData(r).E = [diff(cObsData(r).E)]/Dt;
    %ObsData(r).Z = [diff(cObsData(r).Z)]/Dt;
end
%}

%% Rotate data/synth as needed, scale
%get angle from XY to NE
PAuxFiles = '/Users/bcbirkel/Documents/Seismofiles/hercules_preprocessing/MatlabCodes/AuxFiles/';
FileID = fopen([PAuxFiles 'XYtoLonLat.txt']);
Zone = fscanf(FileID, '%f', 1); An = fscanf(FileID, '%f', 1);

% Select some stations
%iS = [1 3 4 6 7 9 11 13:16 21 23 24 28 31 33 34 37:39 48:50 53 59:61 67 68 70 ...
%      74 77 81 86 88:90 92 94 97 101 107 109];
skip=4;
sta_end=80;
iS = 1:skip:sta_end;
iS = 1:size(ObsData,2);

NF = 1; % Number of figure
NS = 1; %Subplot number
%cd Figures
for r = iS
    % Correct Hercules Moment
    SX1 = MoFactor*(SynS1(r).X);
    SY1 = MoFactor*(SynS1(r).Y);
    SZ1 = MoFactor*(SynS1(r).Z); %TYPICALLY negative, changed for CE
    
    %SX2 = Factor*(SynS2(r).X);
    %SY2 = Factor*(SynS2(r).Y);
    %SZ2 = -Factor*(SynS2(r).Z);
    
    % Rotate X,Y to R,T (skipped)
    %AnR = atan2(Eqk.Stn(r).X-Eqk.X,Eqk.Stn(r).Y-Eqk.Y)*180/pi;
    
    %Syn1(1,:) = cosd(AnR)*SY1+sind(AnR)*SX1;
    %Syn1(2,:) = -sind(AnR)*SY1+cosd(AnR)*SX1;
    %Syn1(3,:) = SZ1;
    
    %Syn2(1,:) = cosd(AnR)*SY2+sind(AnR)*SX2;
    %Syn2(2,:) = -sind(AnR)*SY2+cosd(AnR)*SX2;
    %Syn2(3,:) = SZ2; 
    
    %Rotate X,Y to N,E (where An comes from XYtoLonLat)

    %Syn1(1,:) = cosd(90-An)*SX1+sind(90-An)*SX1;
    %Syn1(2,:) = -sind(90-An)*SX1+cosd(90-An)*SY1;
    Syn1(1,:) = SX1;
    Syn1(2,:) = SY1;
    Syn1(3,:) = SZ1;
    
    % scale data appropriately 
    OD(1,:) = OD_Factor*ObsData(r).N;
    OD(2,:) = OD_Factor*ObsData(r).E;
    OD(3,:) = OD_Factor*ObsData(r).Z;
    
    % Differentiate for velocity
    v_Syn1(1,:) = [diff(Syn1(1,:))];%/DeltaT;
    v_Syn1(2,:) = [diff(Syn1(2,:))];%/DeltaT;
    v_Syn1(3,:) = [diff(Syn1(3,:))];%/DeltaT;
    
    %% Taper, Bandpass and Plot
    t25 = tukeywin(length(v_Syn1(1,:)),0.25);
    %od_t25 = -60:Dt*10:300;
    od_t25 = tukeywin(length(OD(1,:)),0.25); 
    CMP = ['-N';'-E';'-Z'];
    [Am,~] = max(abs(v_Syn1(:)));
    Factor = 1/Am;
    [Am,~] = max(abs(OD(:)));
    ODFactor = 10/Am;
    %ODFactor = 1;
    %TotDist=sqrt((Eqk.Stn(r).X-Eqk.X)^2+(Eqk.Stn(r).Y-Eqk.Y)^2); 
    arclen=distance(Eqk.Stn(r).Lat,Eqk.Stn(r).Lon,Eqk.Lat,Eqk.Lon)
    TotDist=deg2km(arclen)

    starttime = SacHead(r).hd.b;
    starttime = -60;
    timeinc = SacHead(r).hd.delta;
    timeinc = Dt;
    endtime = SacHead(r).hd.npts*SacHead(r).hd.delta;
    endtime = 300;
    OD_T = starttime:timeinc:endtime;
    
    for c = 1:3
        %S1 = bp_bu_co(t25'.*v_Syn1(c,:)*Factor,FMin,FMax,1/Dt,4,2);
        %S2 = bp_bu_co(t25'.*Syn2(c,:),FMin,FMax,1/Dt,4,2); 
        %obd = bp_bu_co(od_t25'.*OD(c,:),FMin,FMax,1/timeinc,5,2);
        %obd = bandpass(od_t25'.*OD(c,:),[1/20,1],1/timeinc);
        %obd = bandpass(OD(c,:),[1/20,1],1/timeinc);
        %test;
        S1 = t25'.*v_Syn1(c,:)*Factor;
        %obd = od_t25'.*OD(c,:);
        obd = OD(c,:);
        
        figure(c);
        %plot(T(1:size(S1,2)),(S1+TotDist),'color',[0 0 0],'linewidth',.3), hold on
        OD_T_trunc=OD_T(1:size(obd,2)); %truncate time to match data length
        plot(OD_T_trunc,((obd.*ODFactor)+SacHead(r).hd.dist),'color',[0.635 0.078 0.184],'linewidth',.3), hold on
        xlabel 'Time (s)', ylabel('Dist from Source (km)'),
        %legend('Synthetic','Observed'), 
        xlim([-60 300]), 
        %ylim(Am*[-1 1])
        %ylim([0 5])
        %text(160,TotDist,[Eqk.Stn(r).Name  CMP(c,:)],'Fontsize', 6), hold on
        text(100,SacHead(r).hd.dist,[SacHead(r).name  CMP(c,:)],'Fontsize', 6, 'Color', 'r'), hold on
    	%if NF>(sta_end/skip-1)
        %    saving = 1
        %    saveas(gcf, ['SeismoDist4' CMP(c,:) '.png'])
            %print(['SeismogramDist.png'],'-bestfit')
        %end
    end
    NF=NF+1;
end


%{
%% UserModification: Path to the Resources directory (catalogs, stations...)
PathRes='/home/scec-00/birkel/hercules_dependencies/Resources/';

%% UserModification: Path to the AuxFiles directory. Here we save or read
%  things that previous codes saved.
PAuxFiles = '/home/scec-00/birkel/hercules_preprocessing/MatlabCodes/AuxFiles/';


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
PathRes='/home/scec-00/birkel/hercules_dependencies/Resources/';

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


% Plot Eqk and stations
for r = iS
    m_plot(Stn(r).Lon,Stn(r).Lat,'v',...
        'markerfacecolor',[0 .45 .74],'markeredgecolor',[0 0 0],...
        'MarkerSize',7,'color','k'); hold on
    m_text(Stn(r).Lon+0.1,Stn(r).Lat,Stn(r).Name)
end


for s=1:size(Eqk,2)
    m_plot(Eqk(s).Loc.Lon,Eqk(s).Loc.Lat,'p','markerfacecolor',...
        [.64 .08 .18],'markeredgecolor',[0 0 0],...
        'MarkerSize',12,'color','k'); hold on
end
%}

