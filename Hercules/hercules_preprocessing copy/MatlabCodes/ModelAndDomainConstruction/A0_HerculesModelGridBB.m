% A0_HerculesModelGrid.m helps to set most of the related model domain 
% files for the Hercules simulation. 
% Written by Alan Juarez (UNAM-USC) Nov, 2016.
% Edited for CyberShake simulation on Nov, 2019.

clear; close all; clc;

% Path to the AuxFiles directory. Here we save things that
% some following codes will need.
PAuxFiles = '/home/scec-00/birkel/hercules_preprocessing/MatlabCodes/AuxFiles/';

% Load mesh coordinate file
FileID = fopen([PAuxFiles 'XYtoLonLat.txt']);

Zone = fscanf(FileID, '%f', 1); Az = fscanf(FileID, '%f', 1);
X0_lat = fscanf(FileID, '%f', 1); Y0_lon = fscanf(FileID, '%f', 1);
XLim(1) = fscanf(FileID, '%f', 1); XLim(2) = fscanf(FileID, '%f', 1);
YLim(1) = fscanf(FileID, '%f', 1); YLim(2) = fscanf(FileID, '%f', 1);
ZLim(1) = fscanf(FileID, '%f', 1); ZLim(2) = fscanf(FileID, '%f', 1);
Corner(1) = fscanf(FileID, '%f', 1); Corner(2) = fscanf(FileID, '%f', 1);
fclose(FileID);

fprintf("check1")

% Computation of the box dimensions, rotation if needed and domain limits
% Y0_lon and X0_lat are in utm coordinates.
X_Lat = XLim(1):1e3:XLim(2);
Y_Lon = YLim(1):1e3:YLim(2);
Z_Depth = [ZLim(1) 100 200 300 400 500 1e3:1e3:ZLim(2)];

fprintf("check1.5")

nP=1;
for i=1:length(X_Lat)
    for j=1:length(Y_Lon)
	fprintf('nP = %f\n', nP)
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

fprintf("check1.75")

% Save XYZ
fileID = fopen([PAuxFiles 'VelModelCoords.xyz'],'w');
for nP=1:length(CoorXY(:,3))
    fprintf(fileID,'%6.2f\t',[CoorXY(nP,1) CoorXY(nP,2) CoorXY(nP,3)]);
    fprintf(fileID,'\n');
end
fclose(fileID);
fprintf("check2")

% Save LonLatZ
fileID = fopen([PAuxFiles 'VelModelCoords.llz'],'w');
for nP=1:length(CoorXY(:,3))
    fprintf(fileID,'%6.2f\t',[CoorLL(nP,2) CoorLL(nP,1) CoorLL(nP,3)]);
    fprintf(fileID,'\n');
end
fclose(fileID);

% Copy files to UCVM and get output
system(['cp /home/scec-00/birkel/hercules_preprocessing/MatlabCodes/AuxFiles/VelModelCoords.llz /home/scec-00/birkel/ucvm-19.4.0/']);
cd /home/scec-00/birkel/ucvm-19.4.0/
system(['./bin/ucvm_query -f ./conf/ucvm.conf -m cvms5 < VelModelCoords.llz > ModelOut.out']);
system(['cp /home/scec-00/birkel/ucvm-19.4.0/ModelOut.out /home/scec-00/birkel/hercules_preprocessing/MatlabCodes/AuxFiles/']);
cd /home/scec-00/birkel/hercules_preprocessing/MatlabCodes/AuxFiles/

fprintf("A0 finished")

% Ploting the grid in LonLatZ and XYZ Coordinates, this is just to show
% more or less the point distribution. Remerber that plots in MatLab are XLon yLat zElev.

%figure;
%stem3(CoorLL(:,2),CoorLL(:,1),CoorLL(:,3)/100000,...
%    '.','linestyle','none','color',rand(1,3))
%xlabel('Lon'), ylabel('Lat'), zlabel('1e5 x Depth [m]')
%set (gca,'Zdir','reverse'), view(2), axis equal

%figure;
%    stem3(CoorXY(:,2)/1000,CoorXY(:,1)/1000,...
%        CoorXY(:,3)/1000,'.','linestyle','none','color',rand(1,3))
%xlabel('Y-Lon [km]'), ylabel('X-Lat [km]'), zlabel('Profundidad [km]')
%set (gca,'Zdir','reverse'), view(2), axis equal

