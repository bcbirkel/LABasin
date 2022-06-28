% Full Waveform Inversion - Code No.C1.
% C1_SourceInput.m creates all the directories and souce related files
% needed for the earthquake forward simulation on Hercules.
% It is very important that in the "Resources Directory" you have a generic
% input.in file with the database and input directories set down and
% with the simulation parameters defined. The only thing that will change
% is the Source Number directory, in the input file use Source999.
% You may also change the input, output and database directories.
% User needs to go througth code and set all the parameters specified with
% "UserModification" legend on the comments.
% Special care must be taken in the folder structure where catalogs,
% resources and databases are located or where they are going to be saved.
% This is probably a mess, but it's helpfull to keep everything in order.

% Written by Alan Juarez (UNAM-USC) Nov, 2016.

clear; close all; clc;

%% UserModification: Path to the AuxFiles directory. Here we save or read
%  things that previous codes saved, for example the database.
PathAux = '/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/hercules_preprocessing/MatlabCodes/AuxFiles/';

% Loading information of the box simulation
% Rotation data. We need to know the number of units
FileID = fopen([PathAux 'XYtoLonLat.txt']);

Zone = fscanf(FileID, '%f', 1); Az = fscanf(FileID, '%f', 1);
fclose(FileID);

% Loading Earthquakes and Stations for the Simulation
load([PathAux 'EarthquakeDB.mat']);

%% UserModification: Reading generic input file and planeinputcoords file
Finp = fileread([PathAux 'inputHPC.in']);

%% UserModification: Reading generic planeinputcoords file
Planes = get_list_files(PathAux, 'planeinputcoords.*');

%% UserModification: Path to Simulation directory
%  Here we are going to create input directories and files
PathSim = '/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/Simulation/Forward/';
mkdir(PathSim);

%% No User Modification Needed: Creating Forward Simulation Files
TD = pwd;

% Slip Function (we assume a point source)
Dt=0.05;
Time=0:Dt:100;

iS=1
cd (PathSim)

% Create Source Directory
namefoldersis = 'Source1';
foldersismo= namefoldersis;
mkdir(foldersismo);
cd(foldersismo);

% OUTPUT
mkdir('output');
cd('output/');
mkdir('forces');
mkdir('planes');
mkdir('stations');
cd ..

% More source parameters
Rake=Eqk(iS).Src.Rake(1);
Strike=Eqk(iS).Src.Str(1)-Az;
Dip=Eqk(iS).Src.Dip(1);
SlipRate=triangularPulse(0,((Eqk(iS).Src.Mo*10^(Eqk(iS).Src.Exp))/(2.2e16))^0.333,Time);
SlipFun=cumtrapz(Time,SlipRate);

% INPUT
mkdir('input');
cd('input/');

fileID = fopen('area.in','w');
fprintf(fileID,'%6.2f\t',1.0);
fclose(fileID);
clear Area fileID;

fileID = fopen('coords.in','w');    % We add datum from model database?
fprintf(fileID,'%6.2f\t',[Eqk(iS).Loc.X Eqk(iS).Loc.Y Eqk(iS).Loc.Depth]);
fclose(fileID);
clear Coords fileID;

fileID = fopen('rake.in','w');
fprintf(fileID,'%i\t',Rake);
fclose(fileID);
clear Rake fileID;

fileID = fopen('strike.in','w');
fprintf(fileID,'%i\t',Strike);
fclose(fileID);
clear Strike fileID;

fileID = fopen('dip.in','w');
fprintf(fileID,'%i\t',Dip);
fclose(fileID);
clear Dip fileID;

fileID = fopen('slip.in','w');
fprintf(fileID,'%6.2f\t',1e20);
fclose(fileID);
clear Slip fileID;

fileID = fopen('slipfunction.in','w');
fprintf(fileID,'%i\t',length(SlipFun));
fprintf(fileID,'%i\t',0);
fprintf(fileID,'%6.4f\t',Dt);
fprintf(fileID,'%6.4f\t',SlipFun/max(SlipFun));
fclose(fileID);
clear fileID;

%Modifying input.in number_output_planes        = 4
S = strrep(Finp,'Forward/Source999',['Forward/' foldersismo]);              % MODIFICAR PATH
S = strrep(S,'number_output_planes        = 4',...
    ['number_output_planes        = ' num2str(size(Planes,1)+1)]); 
fileID = fopen('input.in','w');
fprintf(fileID,S,'*write');
fclose(fileID);
clear S fileID;


%Mesh of the domain
for iP=1:size(Planes,1)
    system(['cp ' PathAux Planes{iP} ' .'])
end

%Set up the  stationcoords.in files. Elevation here is fixed since
%we are going to save them in a plane, given in depth to the topo.
nE=size(Eqk(iS).Stn,2);
fileID = fopen('planeinputcoords.0','w');
fprintf(fileID,'%i\t',[0 1 nE]);
for iA=1:nE
    fprintf(fileID,'\n');
    fprintf(fileID,'%6.2f\t',[Eqk(iS).Stn(iA).X Eqk(iS).Stn(iA).Y 100]);
end

% %% !!!!!!!!!!!!!Requesting an array!!!!!!!!!!!!
% for i=1:11
%     for j=1:11
%         fprintf(fileID,'\n');
%         fprintf(fileID,'%6.2f\t',[700+(i-1)*5 400+(j-1)*5 100]);
%     end
% end
% %% 
fclose(fileID);
clear nE fileID;

cd(TD);
% If you get here without problems then is time to run Hercules for all the
% Sources in the /Simulation/Forward/ directory
