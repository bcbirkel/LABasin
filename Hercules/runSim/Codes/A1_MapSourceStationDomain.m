% Full Waveform Inversion - Code No.2.
% A1_MapSourceStationDomain.m helps to look the station and sources
% catalog and find all of them that are inside the simulation domain. Also
% makes a fancy map with the clasificacion of stations, if the networks are
% especified in the file. This code will give you a lot of problems if you
% want to see a map since mapping needs a lot of specific inputs.
% The most important thing of this code is that it starts creating the
% simulation database with observed data.
% User needs to go througth code and set all the parameters specified with
% "UserModification" legend on the comments.
% Special care must be taken in the folder structure where catalogs,
% resources and databases are located or where they are going to be saved.

clear; close all; clc;

%% UserModification: Path to the Resources directory (catalogs, stations...)
PathRes='/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/Resources/';

%% UserModification: Path to the AuxFiles directory. Here we save or read
%  things that previous codes saved.
PAuxFiles = '/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/hercules_preprocessing/MatlabCodes/AuxFiles/';


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
ColorMap=[loadcmap('tpglpom.c3g'); loadcmap('GrayWhite.c3g')];
m_proj('Miller','lon',[min(BoxSim(:,2))-0.5 max(BoxSim(:,2))+0.5],...
    'lat',[min(BoxSim(:,1))-0.5 max(BoxSim(:,1))+0.5]);
hold on;

% Plot topography and add some other nice things
% m_surf(Lon,Lat,bathymetry,bathymetry),
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


%% User modification needed: Add station network to the database
% Sometimes I see these codes and I feel like a fool
% for complicating my life like this. You can erase all you dont need
FileID = fopen([PathRes '/Catalogs/StationList.txt']);
StN = 1; 
while ~feof(FileID)
    Ntwk = fscanf(FileID, '%s', 1);
    Name = fscanf(FileID, '%s', 1);
    Chan = fscanf(FileID, '%s', 1);
    Lats = fscanf(FileID, '%f', 1);
    Lons = fscanf(FileID, '%f', 1);
    Elev = fscanf(FileID, '%f', 1);
    if StN>1 && length(Chan)>2 && Chan(1)=='B'
        if strcmp(Stn(StN-1).Name,Name) 
            if Chan(3)=='N'
                Stn(StN-1).Cmp(1).Name = Chan;
            elseif Chan(3)=='E'
                Stn(StN-1).Cmp(2).Name = Chan;
            elseif Chan(3)=='Z'
                Stn(StN-1).Cmp(3).Name = Chan;
            end
        else
            Stn(StN).Ntwk = Ntwk;
            Stn(StN).Name = Name;
            Stn(StN).Lon = Lons;
            Stn(StN).Lat = Lats;
            Stn(StN).Elev = Elev;
            if Chan(3)=='N'
                Stn(StN).Cmp(1).Name = Chan;
            elseif Chan(3)=='E'
                Stn(StN).Cmp(2).Name = Chan;
            elseif Chan(3)=='Z'
                Stn(StN).Cmp(3).Name = Chan;
            end
            StN = StN + 1;
        end
    elseif StN==1 && length(Chan)>2 && Chan(1)=='B'
            Stn(StN).Ntwk = Ntwk;
            Stn(StN).Name = Name;
            Stn(StN).Lon = Lons;
            Stn(StN).Lat = Lats;
            Stn(StN).Elev = Elev;
            if Chan(3)=='N'
                Stn(StN).Cmp(1).Name = Chan;
            elseif Chan(3)=='E'
                Stn(StN).Cmp(2).Name = Chan;
            elseif Chan(3)=='Z'
                Stn(StN).Cmp(3).Name = Chan;
            end
            StN = StN + 1;
    end
    clear Name Chan Lats Lons Elev;
end
fclose(FileID);

%% User modification needed: Change directory and CMT catalog file name.
% CMT_Catalog contains all the events you may want to consider
% following the next file format:
% Columns: lon lat str1 dip1 rake1 str2 dip2 rake2 sc iexp name
% rows for all events. Looking in the "GlobalCMT_Project" website choose the
% "GMT psvelomeca input" output format to make it easier.
CMT = fopen([PathRes 'Catalogs/CMTCatalog.txt']);


%% No Modifications needed in this section.
% We are going to read all sources and then delate those that are not
% inside our simulation box
SrN = 1;
Eqk(SrN).Loc.Lon = fscanf(CMT, '%f', 1);
while ~isempty(Eqk(SrN).Loc.Lon)
    Eqk(SrN).Loc.Lat = fscanf(CMT, '%f', 1);
    Eqk(SrN).Loc.Depth = fscanf(CMT, '%f', 1);
    Eqk(SrN).Loc.Depth = 1000*Eqk(SrN).Loc.Depth; %Meters
    Eqk(SrN).Src.Str(1)=fscanf(CMT, '%f', 1);
    Eqk(SrN).Src.Dip(1)=fscanf(CMT, '%f', 1);
    Eqk(SrN).Src.Rake(1)=fscanf(CMT, '%f', 1);
    Eqk(SrN).Src.Str(2)=fscanf(CMT, '%f', 1);
    Eqk(SrN).Src.Dip(2)=fscanf(CMT, '%f', 1);
    Eqk(SrN).Src.Rake(2)=fscanf(CMT, '%f', 1);
    Eqk(SrN).Src.Mo=fscanf(CMT, '%f', 1);
    Eqk(SrN).Src.Exp=fscanf(CMT, '%f', 1);
    Eqk(SrN).Src.Exp=Eqk(SrN).Src.Exp-7;
    Eqk(SrN).Loc.Year=fscanf(CMT, '%f', 1);
    Eqk(SrN).Loc.Month=fscanf(CMT, '%f', 1);
    Eqk(SrN).Loc.Day=fscanf(CMT, '%f', 1);
    Eqk(SrN).Loc.Hour=fscanf(CMT, '%f', 1);
    Eqk(SrN).Loc.Min=fscanf(CMT, '%f', 1);
    Eqk(SrN).Loc.Sec=fscanf(CMT, '%f', 1);
    SrN = SrN + 1;
    Eqk(SrN).Loc.Lon = fscanf(CMT, '%f', 1);
end
Eqk(SrN) = [];
fclose(CMT);

%% No Modifications needed in this section.
% Here we are going to pick all sources and stations inside the simulation
% box. It is probably that something wont look fine when plotting map
iSt=1; CA=[0 .45 .74]; 
for StN=1:size(Stn,2)
    [Y,X] = ll2utm(Stn(StN).Lat,Stn(StN).Lon,Zone);
    CoorXY=[[cosd(Az) -sind(Az); sind(Az) cosd(Az)]*[Y-Y0_lon; X-X0_lat]]';
    if max(XLim)>CoorXY(2)+1e3 && CoorXY(2)-1e3>min(XLim)  ...
            && max(YLim)>CoorXY(1)+1e3 && CoorXY(1)>min(YLim)-1e3
            m_plot(Stn(StN).Lon,Stn(StN).Lat,'v',...
                'markerfacecolor',CA,'markeredgecolor',[0 0 0],...
                'MarkerSize',7,'color','k'); hold on
        StnUsf(iSt)=StN;
        iSt=iSt+1;
    end
    clear X Y CoorXY
end

iSr=1; CS=[.64 .08 .18];
for SrN=1:size(Eqk,2)
    [Y,X] = ll2utm(Eqk(SrN).Loc.Lat,Eqk(SrN).Loc.Lon,Zone);
    CoorXY=[[cosd(Az) -sind(Az); sind(Az) cosd(Az)]*[Y-Y0_lon; X-X0_lat]]';
    if max(XLim)>CoorXY(2)+1e3 && CoorXY(2)-1e3>min(XLim)  ...
            && max(YLim)>CoorXY(1)+1e3 && CoorXY(1)-1e3>min(YLim)
        m_plot(Eqk(SrN).Loc.Lon,Eqk(SrN).Loc.Lat,'p','markerfacecolor',...
            CS,'markeredgecolor',[0 0 0],...
            'MarkerSize',12,'color','k'); hold on
        SrcUsf(iSr)=SrN;
        iSr=iSr+1;
    end
    clear X Y CoorXY ;
end


%% No Modifications needed in this section.
% Plot legends and some other fancy things to the MAP
MinLat = min(BoxSim(:,1));
MinLon = min(BoxSim(:,2));
m_plot(MinLon+0.2,MinLat+0.4,'p',...
    'markerfacecolor',CS,'markeredgecolor',[0 0 0],...
    'MarkerSize',10,'color','k'); hold on
m_plot(MinLon+0.2,MinLat+0.2,'v',...
    'markerfacecolor',CA,'markeredgecolor',[0 0 0],...
    'MarkerSize',8,'color','k'); hold on

m_text(MinLon+0.2,MinLat+0.4,'  Epicenters','FontSize',10)
m_text(MinLon+0.2,MinLat+0.2,'  Receivers','FontSize',10)

%% No Modifications needed in this section.
% Saving station and earthquakes database
Data.Stn=Stn(StnUsf);
Data.Eqk=Eqk(SrcUsf);

save([PAuxFiles 'EqkStnDB.mat'],'Data');


