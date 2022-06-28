% Full Waveform Inversion - Code No.1. 
% A0_HerculesModelGrid.m help to set most of the related model domain files
% for the Hercules simulation. 
% Written by Alan Juarez (UNAM-USC) Nov, 2016.
% Mod Brianna Birkel Sept 2021

clear; close all; clc;

% Path to the AuxFiles directory. Here we save things that
% some following codes will need.
PAuxFiles = '/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/hercules_preprocessing/MatlabCodes/AuxFiles/';
FileID = fopen([PAuxFiles 'XYtoLonLat.txt']);

Zone = fscanf(FileID, '%f', 1); Az = fscanf(FileID, '%f', 1);
X0_lat = fscanf(FileID, '%f', 1); Y0_lon = fscanf(FileID, '%f', 1);
XLim(1) = fscanf(FileID, '%f', 1); XLim(2) = fscanf(FileID, '%f', 1);
YLim(1) = fscanf(FileID, '%f', 1); YLim(2) = fscanf(FileID, '%f', 1);
ZLim(1) = fscanf(FileID, '%f', 1); ZLim(2) = fscanf(FileID, '%f', 1);
Corner(1) = fscanf(FileID, '%f', 1); Corner(2) = fscanf(FileID, '%f', 1);
fclose(FileID);


% Computation of the box dimensions, rotation if needed and domain limits
% Y0_lon and X0_lat are in utm coordinates.
X_Lat=0e3:5e3:760e3;
Y_Lon=0e3:5e3:475e3;
Z_Depth=1e3:4e3:49e3;

nP=1;
for i=1:length(X_Lat)
    for j=1:length(Y_Lon)
        ROT=[cosd(Az) sind(Az); -sind(Az) cosd(Az)]*[Y_Lon(j); X_Lat(i)]+[Y0_lon; X0_lat];
        CLatLon = utm2ll(ROT(1),ROT(2),Zone);
        for k = 1:length(Z_Depth)
            CoorXY(nP,:) = [X_Lat(i) Y_Lon(j) Z_Depth(k)];
            CoorLL(nP,:) = [CLatLon(1) CLatLon(2) Z_Depth(k)]; % Coords in Lat Lot Depth
            nP=1+nP;
        end
        clear ROT CLatLon;
    end
end

  
% Ploting the grid in LonLatZ and XYZ Coordinates, this is just to show
% more or less the point distribution. Remerber that plots in MatLab are XLon yLat zElev.
figure;
stem3(CoorLL(:,2),CoorLL(:,1),CoorLL(:,3)/100000,...
    '.','linestyle','none','color',rand(1,3))
xlabel('Lon'), ylabel('Lat'), zlabel('1e5 x Depth [m]')
set (gca,'Zdir','reverse'), view(2), 

figure;
    stem3(CoorXY(:,2)/1000,CoorXY(:,1)/1000,...
        CoorXY(:,3)/1000,'.','linestyle','none','color',rand(1,3))
xlabel('Y-Lon [km]'), ylabel('X-Lat [km]'), zlabel('Profundidad [km]')
set (gca,'Zdir','reverse'), view(2),

fileID = fopen([PAuxFiles 'GridCoords.xyz'],'w');
for nP=1:length(CoorXY(:,3))
    fprintf(fileID,'%6.4f\t',[CoorLL(nP,1) CoorLL(nP,2) CoorLL(nP,3)]);
    fprintf(fileID,'\n');
end
fclose(fileID);

fileID = fopen([PAuxFiles 'ModelCoords.in'],'w');
fprintf(fileID,'%6.0f',length(CoorXY(:,3)));
fprintf(fileID,'\n');
for nP=1:length(CoorXY(:,3))
    fprintf(fileID,'%6.4f\t',[CoorXY(nP,2) CoorXY(nP,1) CoorXY(nP,3)]);
    fprintf(fileID,'\n');
end
fclose(fileID);

% Saving auxiliary files, planeinputcoords contains the points of the grid
% where deformations will be saved by Hercules.
nPlane=ceil(size(CoorXY,1)/1e5);
lPl=ceil(size(CoorXY,1)/nPlane);

p=1;
for n=1:nPlane
    for m=1:3
        fileID = fopen([PAuxFiles 'planeinputcoords.' num2str(p)],'w');
        fprintf(fileID,'%i\t',m);
        fprintf(fileID,'%i\t',1);
        if n*lPl <= size(CoorXY,1)
            fprintf(fileID,'%i\t',lPl);
            for nP=(n-1)*lPl+1:n*lPl
                fprintf(fileID,'\n');
                fprintf(fileID,'%6.4f\t',[CoorXY(nP,1) CoorXY(nP,2) CoorXY(nP,3)]);
            end
        else
            fprintf(fileID,'%i\t',size(CoorXY,1)-(n-1)*lPl);
            for nP=(n-1)*lPl+1:size(CoorXY,1)
                fprintf(fileID,'\n');
                fprintf(fileID,'%6.4f\t',[CoorXY(nP,1) CoorXY(nP,2) CoorXY(nP,3)]);
            end
        end
        p=p+1;
    end
end
fclose(fileID);

