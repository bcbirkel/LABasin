% C1_SourceInput.m creates all the directories and souce related files
% needed for the earthquake forward simulation on Hercules.
% It is very important that you have a generic
% input.in file with the database and input directories set down and
% with the simulation parameters defined. The only thing that will change
% is the Source Number directory, in the input file use Source999.
% You may also change the input, output and database directories.
% Written by Alan Juarez (UNAM-USC) Nov, 2016.

clear; close all; clc;

% Loading information of the simulation box
FileID = fopen('/home/scec-00/alanjuar/EarthquakePhysics/Codes/AuxFiles/XYtoLonLat.txt');

Zone = fscanf(FileID, '%f', 1); Az = fscanf(FileID, '%f', 1);
X0_lat = fscanf(FileID, '%f', 1); Y0_lon = fscanf(FileID, '%f', 1);
XLim(1) = fscanf(FileID, '%f', 1); XLim(2) = fscanf(FileID, '%f', 1);
YLim(1) = fscanf(FileID, '%f', 1); YLim(2) = fscanf(FileID, '%f', 1);
ZLim(1) = fscanf(FileID, '%f', 1); ZLim(2) = fscanf(FileID, '%f', 1);
Corner(1) = fscanf(FileID, '%f', 1); Corner(2) = fscanf(FileID, '%f', 1);
fclose(FileID);

% Reading generic input file and planeinputcoords file
Finp = fileread('/home/scec-00/alanjuar/EarthquakePhysics/Codes/AuxFiles/inputHPC.in');

% Generate planeinputcoords file
% Saving auxiliary files, planeinputcoords contains the points of the grid
% where deformations will be saved by Hercules.
X_Lat = XLim(1):2e3:XLim(2);
Y_Lon = YLim(1):2e3:YLim(2);
Z_Depth = 50;

CoorXY = zeros(length(X_Lat)*length(Y_Lon)*length(Z_Depth),3);
nP=1;
for i=1:length(X_Lat)
    for j=1:length(Y_Lon)
        for k = 1:length(Z_Depth)
            CoorXY(nP,:) = [X_Lat(i) Y_Lon(j) Z_Depth(k)];
            nP=1+nP;
        end
    end
end

nPlane=ceil(size(CoorXY,1)/1e5);
lPl=ceil(size(CoorXY,1)/nPlane);

p = 0;
for n = 1:nPlane
    fileID = fopen(['/home/scec-00/alanjuar/EarthquakePhysics/Codes/AuxFiles/planeinputcoords.' num2str(p)],'w');
    fprintf(fileID,'%i\t',0);
    fprintf(fileID,'%i\t',1);
    if n*lPl <= size(CoorXY,1)
        fprintf(fileID,'%i\t',lPl);
        for nP = (n-1)*lPl+1:n*lPl
            fprintf(fileID,'\n');
            fprintf(fileID,'%6.2f\t',[CoorXY(nP,1) CoorXY(nP,2) CoorXY(nP,3)]);
        end
    else
        fprintf(fileID,'%i\t',size(CoorXY,1)-(n-1)*lPl);
        for nP = (n-1)*lPl+1:size(CoorXY,1)
            fprintf(fileID,'\n');
            fprintf(fileID,'%6.2f\t',[CoorXY(nP,1) CoorXY(nP,2) CoorXY(nP,3)]);
        end
    end
    p=p+1; 
    fclose all;
end


% Path to Simulation directory
%  Here we are going to create input directories and files
PathSim = '/home/scec-00/alanjuar/EarthquakePhysics/Simulations/';

% Creating Simulation Files
TD = pwd;

% The Dt and Time here are the same as in the input.in file
DtSim = 0.025;
TimeSim = 0:DtSim:200;

% Create random sources
Strike = 360*rand(10,1);
Dip = 90*rand(10,1);
Rake = 360*rand(10,1)-180;

SlipRate = gampdf(TimeSim,2,1.5);
SlipFun = cumtrapz(TimeSim,SlipRate);

M0 = 3.8301e+16; 
X = 200e3*rand(10,1);
Y = 100e3*rand(10,1); 
Z = 20e3*rand(10,1);

%% Create Sources
for ns = 1:10
    % Go to Simulations folder
    cd (PathSim)
    
    % Create Source Directory
    foldersismo = ['Source2_' num2str(ns)];
    mkdir(foldersismo);
    cd(foldersismo);
    
    % create output directories
    mkdir('output');
    cd('output/');
    mkdir('forces');
    mkdir('planes');
    mkdir('stations');
    cd ..
    
    % INPUT
    mkdir('input');
    cd('input/');
    
    fileID = fopen('area.in','w');
    fprintf(fileID,'%6.2f\t',1.0);
    fclose(fileID);
    
    fileID = fopen('coords.in','w');    % We add datum from model database?
    fprintf(fileID,'%6.2f\t',[X(ns) Y(ns) Z(ns)]);
    fclose(fileID);
    
    fileID = fopen('rake.in','w');
    fprintf(fileID,'%i\t',Rake(ns));
    fclose(fileID);
    
    fileID = fopen('strike.in','w');
    fprintf(fileID,'%i\t',Strike(ns));
    fclose(fileID);
    
    fileID = fopen('dip.in','w');
    fprintf(fileID,'%i\t',Dip(ns));
    fclose(fileID);
    
    fileID = fopen('slip.in','w');
    fprintf(fileID,'%6.2f\t',M0);
    fclose(fileID);
    
    fileID = fopen('slipfunction.in','w');
    fprintf(fileID,'%i\t',length(SlipFun));
    fprintf(fileID,'%i\t',0);
    fprintf(fileID,'%6.4f\t',DtSim);
    fprintf(fileID,'%6.4f\t',SlipFun);
    fclose(fileID);
    
    %Modifying input.in number_output_planes        = 4
    S = strrep(Finp,'Simulations/Source999',['Simulations/' foldersismo]);              % MODIFICAR PATH
    S = strrep(S,'number_output_planes        = 4',...
        ['number_output_planes        = ' num2str(nPlane)]);
    fileID = fopen('input.in','w');
    fprintf(fileID,S,'*write');
    fclose(fileID);
    clear S fileID;
    
    
    %Mesh of the domain
    for iP = 1:nPlane
        system(['cp /home/scec-00/alanjuar/EarthquakePhysics/Codes/AuxFiles/planeinputcoords.* .']);
    end
    
    
    cd(TD);
end


% If you get here without problems then is time to run Hercules for all the
% Sources in the /Simulation/ directory
