% Kernel combination
% Written by Alan Juarez (USC) Feb, 2018.
clear; close all; clc;

%% UserModification: Path to the directory where Kernels were saved.
% I recomend a directory where you can preserve them in case you erase sims.
PathAdj = '/home/scec-00/alanjuar/S20100707/Simulation/Adjoint/';

%% UserModification: Path to the AuxFiles directory. Here we save or read
%  things that previous codes saved, for example the database.
PathAux = '/home/scec-00/alanjuar/S20100707/Codes/AuxFiles/';

%% No Modification Needed: Starting plots and calculations

% Model parameters for the kernel
Model = load([PathAux '/Model.out']);
Vp = Model(:,4); Vs = Model(:,5); Rho = Model(:,6);
X = Model(:,2); Y = Model(:,1); Z = Model(:,3);
Grid = load([PathAux '/GridCoords.xyz']);
Lat = Grid(:,1); Lon = Grid(:,2); 
clear Model Grid;


% Initial model for Mu (Shear Modulus) and Kv (Volumetric Constant)
Mu=Rho.*(Vs.^2);
Ka=(Vp.^2).*Rho-(4/3)*Mu;  

nP = 1:9;   % Phase number
nF = 3;   % Phase freq

% Load combination factors
Factor = load([PathAux '/InterFactFrec.mat']);

% Load kernels
for i = 1:length(nP)
    for j=  1:length(nF)
        Phase(i,1).Freq(j,1).Kernels = load([PathAdj 'KernelPhase' num2str(nP(i)) 'Freq' num2str(nF(j)) '.mat']);
    end
end

% Define Kf the objective function
Kf = 0*Vs;
Kf(abs(X-10000) < 10000) = 1;
Kf = [Kf; 0*Kf];

G = zeros([length(Kf) length(nP)*length(nF)]);

k = 1; 
for i = 1:length(nP)
    for j=  1:length(nF)
        G1 = conv([diff(Phase(i,1).Freq(j,1).Kernels.Ks); 0]/1e2,ones(3,1)/3,'same');
        G(:,k) = [Factor.C(i,i,j)*Phase(i,1).Freq(j,1).Kernels.Ks...
            ; G1];
        k = k+1;
    end
end

Alpha = (G'*G)\(G'*Kf);

Alpha = Alpha/max(abs(Alpha));


subplot(2,2,4)
for j=  1:length(nF)
    plot(Alpha,'-..','linewidth',1,'markersize',20), 
    hold on
end
axis 'square', xlim([0 length(nP)+1]), ylim([-1.1 1.1])
xlabel('Phase number'), ylabel('\alpha'), %legend('2 Hz','3 Hz','4 Hz')
title 'Combination coefficient'


% Combine kernels
BASE_C = 0*Factor.BASE(i,:);
Ks = 0*Phase(i,1).Freq(j,1).Kernels.Ks;

k = 1;
for i = 1:length(nP)
    for j=  1:length(nF)
            Ks = Ks+Alpha(k)*Factor.C(i,i,j)*Phase(i,1).Freq(j,1).Kernels.Ks;
            BASE_C = BASE_C+Alpha(k)*Factor.BASE(i,:);
            k = k+1;
    end
end

%% Now you can start plots
    % Creating a Map Mesh for the Map, XY
    XX=unique(X);
    YY=unique(Y);
    ZZ=unique(Z);
    
FileID = fopen([PathAux 'XYtoLonLat.txt']);

Zone = fscanf(FileID, '%f', 1); Az = fscanf(FileID, '%f', 1);
X0_lat = fscanf(FileID, '%f', 1); Y0_lon = fscanf(FileID, '%f', 1);
XLim(1) = fscanf(FileID, '%f', 1); XLim(2) = fscanf(FileID, '%f', 1);
YLim(1) = fscanf(FileID, '%f', 1); YLim(2) = fscanf(FileID, '%f', 1);
fclose(FileID);

iS = 1; % Number of Source
iE = 18; % Number of Stn
%Loading Waveform Database
load([PathAux 'EarthquakeDB.mat']);

    
    %% User Modification Needed: Specify parameters for the plots
    P1=[0 220e3];      % XY coordinates of one point in the the profile (Grid Coordinates)
    P2=[760e3 220e3];  % XY coordinates of other point in the the profile (Grid Coordinates)
    ZMap=9e3;          % Depth of the Map (Will be usefull if you check depths in the A0_Code)
    
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
    
    %% Creating a Map Mesh for the Map, XY and LL
    I = find(abs(Z-ZMap) < 1+min(abs(Z-ZMap)));
    Prop = reshape(Ks(I),length(YY),length(XX))';
    LonP = reshape(Lon(I),length(YY),length(XX))';
    LatP = reshape(Lat(I),length(YY),length(XX))';
    
      
    figure;
    subplot(4,4,[2:11])
    ColorMap=flip(loadcmap('BlueWhiteOrangeRed.c3g'));
    m_proj('Miller','long',[-119.5 -116],...
        'lat',[32.5 36]); hold on;
%     m_proj('Miller','long',[min(BoxSim(:,1)) max(BoxSim(:,1))],...
%         'lat',[min(BoxSim(:,2)) max(BoxSim(:,2))]); hold on;
    m_surf(LonP,LatP,Prop,Prop), shading interp, colormap(ColorMap),
    CaxLi=4*std(Ks); caxis(CaxLi*[-1 1]); hold on; 
    m_grid('box','fancy','tickdir','in','fontsize',9,'fontweight','demi');
%     m_plot(BoxSim(:,1),BoxSim(:,2),'color',[0.5 0.5 0.5],'linewidth',2)
 
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
    PathRes='/home/scec-00/alanjuar/S20100707/Resources/';
    
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
            XYP(nP,:)=[norm([X(n) Y(n)]-P1) Z(n)  Ks(n)];
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
    caxis(CaxLi*[-1 1]); hold on, axis ij, %axis equal
    axis(0.001*[hh(1)+200e3 hh(end)-200e3 ZZ(1) ZZ(end)])
    xlabel ('Distance [km]'), ylabel ('Depth [km]')
    stem3(0.001*XYS(1),0.001*XYS(2),10e3,'p',...
        'markerfacecolor',[1 0.84 0],'markeredgecolor',[0 0 0],...
        'MarkerSize',12,'color','k'); hold on
    stem3(0.001*XYR(1),0.001*XYR(2),10e3,'v',...
        'markerfacecolor',[0.47 0.67 0.19],'markeredgecolor',[0 0 0],...
        'MarkerSize',10,'color','k'); hold on
    set (gca,'Xdir','reverse')
    
    set(gcf,'renderer','zbuffer')
    set(gcf,'color','w')
    myaa([8 4])



%% Plot Seismogram Decomposition %%%
dt=0.25;     %Simulations dt
t=0:dt:200;  %Same as simulations

figure;
subplot(2,2,1)
plot(t,sum(Factor.BASE,1)/max(sum(Factor.BASE,1)),'k','linewidth',2), hold on

subplot(2,2,1)
plot(t,BASE_C/max(BASE_C),'r','linewidth',1.5),
legend('Data','Localized','location','southeast')
xlim([0 t(end)]), ylim(1.5*[-1 1]),
title('Localization')
xlabel('Time [s]') , ylabel('Displacement [m]')

