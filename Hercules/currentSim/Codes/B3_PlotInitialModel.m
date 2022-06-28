% Full Waveform Inversion - Code No G1.
% G1_ModelUpdate.m reads the model database and updates the velocities
% using all the computed Kernes stored in the database.
% Special care must be taken in the folder structure where catalogs,
% resources and databases are located or where they are going to be saved.
% This is probably a mess, but it's helpfull to keep everything in order.

% Written by Alan Juarez (UNAM-USC) Dec, 2016.
clear; close all; clc;

%% UserModification: Path to the AuxFiles directory. Here we save or read
%  things that previous codes saved, for example the database.
PathAux = '/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/hercules_preprocessing/MatlabCodes/AuxFiles/';

%% UserModification: Path to the Resources directory (catalogs, stations...)
PathRes='/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/Simulation/Database/';

%Loading Waveform Database
load([PathAux 'EarthquakeDB.mat']);

% Loading information of the box simulation and rotation data.
FileID = fopen([PathAux 'XYtoLonLat.txt']);

Zone = fscanf(FileID, '%f', 1); Az = fscanf(FileID, '%f', 1);
X0_lat = fscanf(FileID, '%f', 1); Y0_lon = fscanf(FileID, '%f', 1);
XLim(1) = fscanf(FileID, '%f', 1); XLim(2) = fscanf(FileID, '%f', 1);
YLim(1) = fscanf(FileID, '%f', 1); YLim(2) = fscanf(FileID, '%f', 1);
fclose(FileID);

iS=1;
iE=18;

FileID = fopen([PathRes 'vel_profiles_unidad1.surf']);
nx=fscanf(FileID, '%i', 1);
ny=fscanf(FileID, '%i', 1);
XX=fscanf(FileID, '%f', nx);
YY=fscanf(FileID, '%f', ny);
fclose(FileID);
for i=1:length(XX)
    for j=1:length(YY)
        ROT=[cosd(Az) sind(Az); -sind(Az) cosd(Az)]*[YY(j); XX(i)]+[Y0_lon; X0_lat];
        [LatP(i,j),LonP(i,j)] = utm2ll(ROT(1),ROT(2),Zone);
        clear ROT 
    end
end   

FileID = fopen([PathRes 'vel_profiles_unidad1.fun']);
NumProf= fscanf(FileID, '%i', 1);
n=1;
for i=1:nx
    for j=1:ny
        nz=fscanf(FileID, '%i', 1);
        X(n:n+nz-1)=XX(i);
        Y(n:n+nz-1)=YY(j);
        Z(n:n+nz-1)=fscanf(FileID, '%f', nz);
        Vp(n:n+nz-1)=fscanf(FileID, '%f', nz);
        Vs(n:n+nz-1)=fscanf(FileID, '%f', nz);
        Rho(n:n+nz-1)=fscanf(FileID, '%f', nz);
        n=n+nz;
    end
end
fclose(FileID);

Vs = Vs/1e3;

%% User Modification Needed: Specify parameters for the plots
P1=[0 220e3];     % XY coordinates of one point in the the profile (Grid Coordinates)
P2=[760e3 220e3];   % XY coordinates of other point in the the profile (Grid Coordinates)
ZMap=1e3;           % Depth of the Map (Will be usefull if you check depths in the A0_Code) 

ZZ=unique(Z);


%% Creating a Map Mesh for the Map, XY and LL
I = find(abs(Z-ZMap) < 1+min(abs(Z-ZMap)));
Prop=reshape(Vs(I),length(YY),length(XX))';

%Map Figure
nP=1;
for i=1:2
    for j=1:2
        ROT=[cosd(Az) sind(Az); -sind(Az) cosd(Az)]*[YLim(i); XLim(j)]+[Y0_lon; X0_lat];
        [BoxSim(nP,2),BoxSim(nP,1)] = utm2ll(ROT(1),ROT(2),Zone);
        nP=nP+1;
    end
end
BoxSim(nP,:)=BoxSim(1,:);
BoxSim([3 4],:)=BoxSim([4 3],:);

figure;
subplot(4,4,[2:11])
ColorMap=flip(loadcmap('BlueWhiteOrangeRed.c3g'));
m_proj('Miller','long',[min(BoxSim(:,1)) max(BoxSim(:,1))],...
    'lat',[min(BoxSim(:,2)) max(BoxSim(:,2))]); hold on;
% m_surf(LonP,LatP,Prop,Prop), shading interp, colormap(ColorMap), 
CaxLi=1.1*min(Vs); CaxLs=min(Vs)+0.9*(max(Vs)-min(Vs));
caxis([CaxLi CaxLs]); hold on; colorbar('west')
m_grid('box','fancy','tickdir','in','fontsize',9,'fontweight','demi'); 
m_plot(BoxSim(:,1),BoxSim(:,2),'color',[0.63 0.07 0.18],'linewidth',2)

% Profile Line
ROT=[cosd(Az) sind(Az); -sind(Az) cosd(Az)]*[P1(2); P1(1)]+[Y0_lon; X0_lat];
[LLP(1,2),LLP(1,1)] = utm2ll(ROT(1),ROT(2),Zone); clear ROT
ROT=[cosd(Az) sind(Az); -sind(Az) cosd(Az)]*[P2(2); P2(1)]+[Y0_lon; X0_lat];
[LLP(2,2),LLP(2,1)] = utm2ll(ROT(1),ROT(2),Zone); clear ROT

m_plot(LLP(:,1),LLP(:,2),'--','linewidth',2,'color',0.3*[1 1 1])

% Plot Estaciones y Fuentes Sismos
for iR1=iE%1:size(Eqk(iS).Stn,2)
    m_plot(Eqk(iS).Stn(iR1).Lon,Eqk(1).Stn(iR1).Lat,'v',...
        'markerfacecolor',[0.47 0.67 0.19],'markeredgecolor',[0 0 0],...
        'MarkerSize',8,'color','k'); hold on
end
m_plot(Eqk(iS).Loc.Lon,Eqk(iS).Loc.Lat,'p',...
    'markerfacecolor',[1 0.84 0],'markeredgecolor',[0 0 0],...
    'MarkerSize',12,'color','k'); hold on

% No modification needed in this secction
% loading topography file
PathRes='/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/Resources/';

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

[Lon,Lat]=meshgrid(lon,lat);    % lon and lat come from elevationdataforaxis.mat
m_contour(lon,lat,bathymetry,[0, 0],'Color',[0.2 0.2 0.2],'linewidth',2);
load([PathRes 'RollinsData/caldata.mat']);
m_plot(uclafaults(:,1),uclafaults(:,2),'color',[0.6 0.6 0.6],'linewidth',1)

%% Creating Profile
% Profile parameters (Fancy Things Designed by Alan) 
AL=-(P2(2)-P1(2));
BL=(P2(1)-P1(1));
CL=(P2(2)-P1(2))*P1(1)-(P2(1)-P1(1))*P1(2);

% Look for the points in the profile
nP=1;
for n=1:length(X)
    if abs(AL*X(n)+BL*Y(n)+CL)/(AL^2+BL^2)^(1/2) < 2.5e3 
        XYP(nP,:)=[norm([X(n) Y(n)]-P1) Z(n)  Vs(n)];
        nP=nP+1;
    end
end

% Mesh for the profile
[hh,zz]=meshgrid([min(XYP(:,1)):min(diff(unique(Z))):max(XYP(:,1))],[min(XYP(:,2)):min(diff(unique(Z))):max(XYP(:,2))]);

% Here we need interpolation since the points are not regular distributed
Prop = griddata(XYP(:,1),XYP(:,2),XYP(:,3),hh,zz,'cubic'); %Interpolamos
XYS=[norm([Eqk(iS).Loc.X Eqk(iS).Loc.Y]-[0 -CL/BL]) Eqk(iS).Loc.Depth];
XYR=[norm([Eqk(iS).Stn(iE).X Eqk(iS).Stn(iE).Y]-[0 -CL/BL]) 1e3];

subplot(4,4,[13.75 15.25])
surf(0.001*hh,0.001*zz,Prop), shading interp, view(2),  colormap(ColorMap),
caxis([CaxLi CaxLs]); hold on, axis ij, %axis equal
axis(0.001*[hh(1) hh(end) ZZ(1) ZZ(end)])
xlabel ('Distance (km)        Vs (km/s)'), ylabel ('Depth (km)')
stem3(0.001*XYS(1),0.001*XYS(2),10e3,'p',...
    'markerfacecolor',[1 0.84 0],'markeredgecolor',[0 0 0],...
    'MarkerSize',12,'color','k'); hold on
stem3(0.001*XYR(1),0.001*XYR(2),10e3,'v',...
    'markerfacecolor',[0.47 0.67 0.19],'markeredgecolor',[0 0 0],...
    'MarkerSize',10,'color','k'); hold on
set (gca,'Xdir','reverse')

set(gcf,'renderer','zbuffer')
set(gcf,'color','w')
myaa([16 2])
print(['./Figures/VelocityModel.pdf'],'-dpdf')
    
    
    
