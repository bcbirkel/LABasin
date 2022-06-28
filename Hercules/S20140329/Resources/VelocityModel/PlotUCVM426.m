% Process and plot CVM426
% Written by Alan Juarez (UNAM-USC) Dec, 2016.
clear; close all; clc;

ZMap=500;

FileID = fopen('cvms426_final_model.txt');
Header=fgetl(FileID);

N = 1;  Depth(N)=ZMap;
while ~feof(FileID) 
    X(N,1)=fscanf(FileID, '%i', 1);
    if ~isempty(X)
        Y(N,1)=fscanf(FileID, '%i', 1);
        Z(N,1)=fscanf(FileID, '%i', 1);
        UTM_X(N,1)=fscanf(FileID, '%f', 1);
        UTM_Y(N,1)=fscanf(FileID, '%f', 1);
        Lon(N,1)=fscanf(FileID, '%f', 1);
        Lat(N,1)=fscanf(FileID, '%f', 1);
        Depth(N,1)=fscanf(FileID, '%f', 1);
        Vp(N,1)=fscanf(FileID, '%f', 1);
        Vs(N,1)=fscanf(FileID, '%f', 1);
        Rho(N,1)=fscanf(FileID, '%f', 1);
        
        if Depth(N)-ZMap > 1
            break;
        end
        if abs(Depth(N)-ZMap)<1
            N = N + 1;
        end
        
    end
end
fclose(FileID);

x=unique(X);
y=unique(Y);

LonP=reshape(Lon(1:end-1),length(x),length(y))';
LatP=reshape(Lat(1:end-1),length(x),length(y))';
Prop=reshape(Vs(1:end-1),length(x),length(y))';


figure;
subplot(4,1,[1 3])
ColorMap=flip(loadcmap('Spectrum120.c3g'));
m_proj('Miller','long',[min(Lon) max(Lon)],...
    'lat',[min(Lat) max(Lat)]); hold on;
m_surf(LonP,LatP,Prop,Prop), shading interp, colormap(ColorMap), 
CaxLi=0.9*min(Vs); CaxLs=1.1*max(Vs);
caxis([CaxLi CaxLs]); hold on; colorbar('west')
m_grid('box','fancy','tickdir','in','fontsize',12,'fontweight','demi'); 

% loading topography file
PathRes='/home/dude/alanjuar/Documents/LaHabra/Resources/';

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

m_contour(lon,lat,bathymetry,[0, 0],'Color',[0.2 0.2 0.2],'linewidth',2);
load([PathRes 'RollinsData/caldata.mat']);
m_plot(uclafaults(:,1),uclafaults(:,2),'color',[0.6 0.6 0.6],'linewidth',1)

