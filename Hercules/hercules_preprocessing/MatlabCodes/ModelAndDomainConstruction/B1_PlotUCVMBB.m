% Process and plot CVM426
% Written by Alan Juarez (UNAM-USC) Dec, 2016.
clear; close all; clc;

% Load XYZ mesh data
Mesh = load(['/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/hercules_preprocessing/MatlabCodes/AuxFiles/VelModelCoords.xyz']);
X = Mesh(:,1);
Y = Mesh(:,2);
clear Mesh;

% Load Earthquake database
load /Users/bcbirkel/Documents/GitHub/LABasin/Hercules/hercules_preprocessing/MatlabCodes/DataPreparation/EarthquakeDatabase.mat;

% Loading information of the simulation box
FileID = fopen('/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/hercules_preprocessing/MatlabCodes/AuxFiles/XYtoLonLat.txt');

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

% Load UCVM model
fileID = fopen('/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/hercules_preprocessing/MatlabCodes/AuxFiles/ModelOut.out','r');
formatSpec = ['%f %f %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f'];
A = textscan(fileID,formatSpec);
fclose('all');

% Assign properties to arrays
Lon = A{1,1};
Lat = A{1,2};
Z = A{1,3};
% Surf = A{1,4};
% vs30 = A{1,5};
% crustal = A{1,6};
% cr_vp = A{1,7};
% cr_vs = A{1,8};
% cr_rho = A{1,9};
% gtl = A{1,10};
% gtl_vp = A{1,11};
% gtl_vs = A{1,12};
% gtl_rho = A{1,13};
% cmb_algo = A{1,14};
Vp = A{1,15};
Vs = A{1,16};
Rho = A{1,17};
clear A;

% % % Estimate 30% velocity drop
% In = find(Vs <= 500 & Z == 0);
% Vs(In) = 0.7*Vs(In);


% Select map surface
ZMap = 0.5e3;
In = find(Z == ZMap);

% unique arrays
x = unique(X);
y = unique(Y);

% Property to plot
Prop = reshape(Vs(In),length(y),length(x))';
Prop = imgaussfilt(0.001*Prop,1);

LonP=reshape(Lon(In),length(y),length(x))';
LatP=reshape(Lat(In),length(y),length(x))';

figure;
m_proj('Miller','lon',[min(BoxSim(:,2))-0.1 max(BoxSim(:,2))+0.1],...
    'lat',[min(BoxSim(:,1))-0.1 max(BoxSim(:,1))+0.1]); hold on;

m_coast('patch',[160 160 160]/255,'edgecolor','k'); 
m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor',[204 226 255]/255);

m_pcolor(LonP,LatP,Prop); caxis([0 4.5]);
cm = colorbar('eastoutside'); title(cm,'Vs (km/s)');
cm.Limits = [0.5 4];
colormap(flip(loadcmap('WhiteBlueGreenYellowRed.c3g')));

% loading topography file (if got it from ETOPO you are done!)
elev=load ('/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/Resources/Topography/etopo1_bedrock.xyz');
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
m_contour(lon,lat,bathymetry,[0, 0],'Color',[0.2 0.2 0.2],'linewidth',2);


% Load and plot faults
load('/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/Resources/RollinsData/caldata.mat');
m_plot(uclafaults(:,1),uclafaults(:,2),'color',[0.4 0.4 0.4],'linewidth',1)
% m_plot(stateborders(:,2),stateborders(:,1),'color',[0.1 0.1 0.1],'linewidth',2)


% Add grid and plot simulation box
m_plot(BoxSim(:,2),BoxSim(:,1),'color',[0.27 0.27 0.27],'linewidth',2)

% % Plot Eqk and stations
for ns = 1:size(Eqk,2)
    m_plot(Eqk(ns).Lon,Eqk(ns).Lat,'p','markerfacecolor',...
        [.64 .08 .18],'markeredgecolor',[0 0 0],...
        'MarkerSize',12,'color','k'); hold on
    for nr = 1:size(Eqk(ns).Stn,2)
        m_plot(Eqk(ns).Stn(nr).Lon,Eqk(ns).Stn(nr).Lat,'v',...
            'markerfacecolor',[0 .45 .74],'markeredgecolor',[0 0 0],...
            'MarkerSize',5,'color','k'); hold on
    end
end

% Plot legends and some other fancy things to the MAP
MinLat = min(BoxSim(:,1));
MinLon = min(BoxSim(:,2));
m_plot(MinLon+0.2,MinLat+0.3,'p',...
    'markerfacecolor',[.64 .08 .18],'markeredgecolor',[0 0 0],...
    'MarkerSize',10,'color','k'); hold on
m_plot(MinLon+0.2,MinLat+0.2,'v',...
    'markerfacecolor',[0 .45 .74],'markeredgecolor',[0 0 0],...
    'MarkerSize',8,'color','k'); hold on

m_text(MinLon+0.2,MinLat+0.3,'  Epicenters','FontSize',10)
m_text(MinLon+0.2,MinLat+0.2,'  Receivers','FontSize',10)

% print -dpng -r800 -painters CVMS5_500m
saveas(gcf,'CVMS5_500mBB.png')
