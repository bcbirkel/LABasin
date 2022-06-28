% A1_MapSourceStationDomain.m helps to look the station and sources
% catalog and find all of them that are inside the simulation domain.
% Special care must be taken in the folder structure where catalogs,
% resources and databases are located or where they are going to be saved.

clear; close all; clc;

% Load Earthquake database
load /home/scec-00/birkel/hercules_preprocessing/MatlabCodes/DataPreparation/EarthquakeDatabase.mat;

% Loading information of the simulation box
FileID = fopen('/home/scec-00/birkel/hercules_preprocessing/MatlabCodes/AuxFiles/XYtoLonLat.txt');

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

% load topography for the map
% In my case I have fancy files with topography, you can edit this
% seccion according to your needs. Topography can be easyly included by 
% downloading files from ETOPO1.

% loading topography file (if got it from ETOPO you are done!)
elev=load ('/home/scec-00/birkel/hercules_dependencies/Resources/Topography/etopo1_bedrock.xyz');
lon = sort(unique(elev(:,1)));
lat = sort(unique(elev(:,2)),'descend');
i=1;
for jj=1:length(lat)
    for ii=1:length(lon)
        bathymetry(jj,ii)=elev(i,3);
        i=i+1;
    end
end

% Map setting
[Lon,Lat]=meshgrid(lon,lat);    
% Creating and setting Map parameters. Here we are using information in the
% resources directory
figure;
ColorMap=[loadcmap('tpglpom.c3g'); loadcmap('GrayWhite.c3g')];
m_proj('Miller','lon',[min(BoxSim(:,2))-0.1 max(BoxSim(:,2))+0.1],...
    'lat',[min(BoxSim(:,1))-0.1 max(BoxSim(:,1))+0.1]);
hold on;

% Plot topography and add some other nice things
m_pcolor(Lon,Lat,bathymetry); 
view (2), colormap(ColorMap),
caxis([-0.5*max(max(bathymetry)) 0.5*max(max(bathymetry))]), shading interp
m_contour(lon,lat,bathymetry,[0, 0],'Color',[0.2 0.2 0.2],'linewidth',1);

% Load and plot faults
load('/home/scec-00/birkel/hercules_dependencies/Resources/RollinsData/caldata.mat');
m_plot(uclafaults(:,1),uclafaults(:,2),'color',[0.5 0.5 0.5],'linewidth',1)
% m_plot(stateborders(:,2),stateborders(:,1),'color',[0.1 0.1 0.1],'linewidth',2)


% Add grid and plot simulation box
m_grid('tickdir','in','box','fancy','fontsize',12,'fontweight','demi');
m_plot(BoxSim(:,2),BoxSim(:,1),'color',[0.27 0.27 0.27],'linewidth',3)


% Here we are going to pick all sources and stations inside the simulation
% box.

% Initiate counters to delate Eqks and Stns
iSr=1; CS=[.64 .08 .18];

for ns = 1:size(Eqk,1)
    
    % Estimate UTM coordinates
    [Y,X] = ll2utm(Eqk(ns).Lat,Eqk(ns).Lon,Zone);
    
    % Rotate to domain
    CoorXY=([cosd(Az) -sind(Az); sind(Az) cosd(Az)]*[Y-Y0_lon; X-X0_lat])';
    
    % Is it within the domain?
    if XLim(2) > CoorXY(2)+1e3 && CoorXY(2)-1e3 > XLim(1)  ...
            && YLim(2) > CoorXY(1)+1e3 && CoorXY(1)-1e3 > YLim(1)
            m_plot(Eqk(ns).Lon,Eqk(ns).Lat,'p','markerfacecolor',...
                CS,'markeredgecolor',[0 0 0],...
                'MarkerSize',12,'color','k'); hold on
            Eqk(ns).X = CoorXY(2);
            Eqk(ns).Y = CoorXY(1);
    else
        SrcUsf(iSr) = ns;
        iSr=iSr+1;
    end
    clear X Y CoorXY ;
    
    % Iterate over stations
    iSt=1; CA=[0 .45 .74];  % Stations to delete counter
    for nr = 1:size(Eqk(ns).Stn,2)
        
        % Estimate UTM coordinates
        [Y,X] = ll2utm(Eqk(ns).Stn(nr).Lat,Eqk(ns).Stn(nr).Lon,Zone);
        
        % Rotate to domain
        CoorXY=([cosd(Az) -sind(Az); sind(Az) cosd(Az)]*[Y-Y0_lon; X-X0_lat])';
        
        % Is it within the domain?
        if XLim(2) > CoorXY(2)+1e3 && CoorXY(2)-1e3 > XLim(1)  ...
            && YLim(2) > CoorXY(1)+1e3 && CoorXY(1)-1e3 > YLim(1)
                m_plot(Eqk(ns).Stn(nr).Lon,Eqk(ns).Stn(nr).Lat,'v',...
                    'markerfacecolor',CA,'markeredgecolor',[0 0 0],...
                    'MarkerSize',7,'color','k'); hold on
                Eqk(ns).Stn(nr).X = CoorXY(2);
                Eqk(ns).Stn(nr).Y = CoorXY(1);
         else
            StnUsf(iSt) = nr;
            iSt=iSt+1;
        end
        clear X Y CoorXY
    end
    Eqk(ns).Stn(StnUsf) = [];
    clear StnUsf;
end
% Eqk(SrcUsf) = [];

% Plot legends and some other fancy things to the MAP
MinLat = min(BoxSim(:,1));
MinLon = min(BoxSim(:,2));
m_plot(MinLon+0.2,MinLat+0.4,'p',...
    'markerfacecolor',CS,'markeredgecolor',[0 0 0],...
    'MarkerSize',10,'color','k'); hold on
m_plot(MinLon+0.2,MinLat+0.2,'v',...
    'markerfacecolor',CA,'markeredgecolor',[0 0 0],...
    'MarkerSize',8,'color','k'); hold on

m_text(MinLon+0.2,MinLat+0.4,'  Epicenter','FontSize',10)
m_text(MinLon+0.2,MinLat+0.2,'  Receivers','FontSize',10)

% Saving station and earthquakes database
save('/home/scec-00/birkel/hercules_preprocessing/MatlabCodes/DataPreparation/EarthquakeDatabase.mat','Eqk');

