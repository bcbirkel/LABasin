% Full Waveform Inversion - Code No.E1-Test
% E0_InputAdjointPhases.m creates all the directories and souce related files
% needed for the earthquake ADJOINT simulation on Hercules.
% It is very important that in the "Resources Directory" you have a generic
% input.in file with the database and input directories setted down and
% with the simulation parameters defined. The only thing that will change
% is the Source Number directory, in the input file use Source999.

% This code should be run after all the Hercules simulations are done.

% This code is to set files for individual kernels. It helps to chek in the
% begining that everything is going well.
% This code calls the GSDF function. Add the path.

% User needs to go througth code and set all the parameters specified with
% "UserModification" legend on the comments.
% Special care must be taken in the folder structure where catalogs,
% resources and databases are located or where they are going to be saved.
% This is probably a mess, but it's helpfull to keep everything in order.

% Written by Alan Juarez (UNAM-USC) Jan, 2016.
clear; close all; clc;

%% UserModification: Path to the Source directory. In this directory are
% saved the outputs from the Hercules Simulation
PathSrc = '/home/dude/alanjuar/S20100707/Simulation/ForwardStn45/';

%% UserModification: Path to the Simulation directory.
%  Here we are going to create input directories and files
PathSim = '/home/dude/alanjuar/S20100707/Simulation/AdjointStn45/';

%% UserModification: Specify the desired minimum and maximum
% frequencies for waveform filtering. (it is a super good idea these match with
% the hercules simulation)
Mo= 4.0654e+20; % Seismic Moment from Hercules
FMax=1/6; %Upper Frequecy
FMin=1/50;  %Lower frequency
dt=0.25; 
t=0:dt:100;
iS=1;   %Source number
nE=45;   %Station number (depends on planeinput.0)
iP=2;   %Component (usually 1=x=lat=r, 2=y=lon=t, 3=z=depth)

%% UserModification: Path to the AuxFiles directory. Here we save or read
%  things that previous codes saved, for example the database.
PathAux = '/home/dude/alanjuar/S20100707/Codes/AuxFiles/';
load([PathAux 'EarthquakeDB.mat']);% Load Earthquakes and Stations
Finp = fileread([PathAux 'inputHPC.in']);

%% No modification needed in this secction
% Loading information of the box simulation
% Rotation data. We are going to set the stations coordinates in the XY
% system.
FileID = fopen([PathAux 'XYtoLonLat.txt']);

fscanf(FileID, '%f', 1); Az = fscanf(FileID, '%f', 1);
fclose(FileID); An=abs(Az);

%% UserModification: Reading generic planeinputcoords file
Planes = get_list_files(PathAux, 'planeinputcoords.*');
TD=pwd;

display(['Source ' num2str(iS)])

%% No modification needed in this secction
% Ploting recording at each station. Open plane
PathPlane=[PathSrc '/Source1/output/planes/planedisplacements.0'];
Plane=fopen(PathPlane);

% Header and auxiliary data
Aux=  fread(Plane, 8,'double');
BoxCorners=reshape(Aux,4,2);
LengthY=  fread(Plane, 1,'double');
LengthX=  fread(Plane, 1,'double');
PointsDip = fread(Plane, 1,'int');
PointsStrike    = fread(Plane, 1,'int');
DeltaT  = fread(Plane, 1,'double');
Steps      = fread(Plane, 1,'int');
TimeSteps = round(t(end)/(Steps*DeltaT)+1);

% Reading Components for all stations
for nM=1:TimeSteps
    for iE=1:PointsStrike
        SynS(iE).X(nM) = fread(Plane, 1,'double');
        SynS(iE).Y(nM) = fread(Plane, 1,'double');
        SynS(iE).Z(nM) = fread(Plane, 1,'double');
    end
end
fclose(Plane);

Factor=(Eqk(iS).Src.Mo*10^Eqk(iS).Src.Exp)/Mo;
SX = Factor*SynS(nE).X;
SY = Factor*SynS(nE).Y;
SZ = Factor*SynS(nE).Z;

DN=Eqk(iS).Stn(nE).Cmpt(1).Signal;
DE=Eqk(iS).Stn(nE).Cmpt(2).Signal;
DZ=Eqk(iS).Stn(nE).Cmpt(3).Signal;

%Take into account if any component has no record
if norm(DZ)==0
    SZ = 0*SZ;
end
if norm(DN)==0 && norm(DE)==0
    SX = 0*SX;
    SY = 0*SY;
end
if norm(DN)==0 || norm(DE)==0
    DXS=SX;
    DYS=SY;
    if norm(DN)==0
        DESR = cosd(An)*DYS-sind(An)*DXS;
        DNSR = 0*(sind(An)*DYS+cosd(An)*DXS);
    elseif norm(DE)==0
        DESR = 0*(cosd(An)*DYS-sind(An)*DXS);
        DNSR = sind(An)*DYS+cosd(An)*DXS;
    end
    SY = cosd(An)*DESR+sind(An)*DNSR;
    SX = -sind(An)*DESR+cosd(An)*DNSR;
    clear DXS DYS DNSR DESR;
end

%Rotate N and E to X and Y
DY=cosd(An)*DE+sind(An)*DN;
DX=-sind(An)*DE+cosd(An)*DN;
clear DN DE;

%Rotate X,Y to R,T
AnR=atan2(Eqk(iS).Stn(nE).X-Eqk(iS).Loc.X,Eqk(iS).Stn(nE).Y-Eqk(iS).Loc.Y)*180/pi;
Obs(1).Signal = cosd(AnR)*DY+sind(AnR)*DX;
Obs(2).Signal = -sind(AnR)*DY+cosd(AnR)*DX;
Obs(3).Signal = DZ;

Syn(1).Signal = cosd(AnR)*SY+sind(AnR)*SX;
Syn(2).Signal = -sind(AnR)*SY+cosd(AnR)*SX;
Syn(3).Signal = SZ;

clear DX DY DZ SX SY SZ;

SynF=bp_bu_co(Syn(iP).Signal,FMin,FMax,1/dt,4,2)';

%%% Make Seismogram Decomposition %%%
subplot(3,2,1)
plot(t,SynF,'k','linewidth',2), hold on

[BASE,CTF] = atom_decm(SynF,dt,.5);
% % Windowing 
% for i=1:t(end)
%     SF = exp(-2^2*(t-t(i/dt+1)-1).^2)';
%     BASE(i,:) = SF.* SynF;
% end

subplot(3,2,1)
plot(t,sum(BASE,1),'r','linewidth',2),
legend('Synthetic','Synthesized','location','southeast')
xlim([0 t(end)]), grid minor, title('Data'), xlabel('Time [s]')
set(gca,'fontsize',12,'fontweight','demi')

MMax=size(BASE,1);
subplot(3,2,[2 4 6])
for m=1:MMax %
    plot(t,m+BASE(m,:)/max(max(0.5*BASE)),'k','linewidth',1.5), hold on
    %     text(1,m+0.12,['U_' num2str(m)],'fontsize',10,'fontweight','demi')
end
xlim([0 t(end)]), ylim([0 m+1]), grid minor,
xlabel 'Time [s]', ylabel 'Atoms'
title('Dictionary'), set(gca,'fontsize',12,'fontweight','demi')



for m=1:MMax 
    %Path to simulation and forlders%
    cd (PathSim)
    NameDir = 'Phase';
    DirSource= strcat(NameDir,num2str(m));
    mkdir(DirSource);
    cd(DirSource);
    
    %%%%%%%%%%% OUTPUT %%%%%%%%%
    mkdir('output');
    cd('output/');
    mkdir('forces');
    mkdir('planes');
    mkdir('stations');
    cd ..
    
    %%%%%%%%%%% INPUT %%%%%%%%%
    mkdir('input');
    cd('input/');
    
    fileID = fopen('weightsxyz.in','w');
    Wrt=[0 0 0];     Wrt(iP)=1;
    
    Wxy(2)=cosd(AnR)*Wrt(1)-sind(AnR)*Wrt(2);
    Wxy(1)=sind(AnR)*Wrt(1)+cosd(AnR)*Wrt(2);
    
    fprintf(fileID,'%i\t',1e20*Wxy);
    fprintf(fileID,'\n');
    fclose(fileID);
    clear fileID;
    
    fileID = fopen('coords.in','w');
    fprintf(fileID,'%6.4f\t',[Eqk(iS).Stn(nE).X Eqk(iS).Stn(nE).Y 100]);
    fprintf(fileID,'\n');
    fclose(fileID);
    clear Coords fileID;
    
    fileID = fopen('sourcefunction.in','w');
    fprintf(fileID,'%i\t',length(BASE(m,:)));
    fprintf(fileID,'%i\t',0);
    fprintf(fileID,'%6.4f\t',dt);
    fprintf(fileID,'%e\t',flip(BASE(m,:)));
    fprintf(fileID,'\n');
    fclose(fileID);
    clear fileID;
    
    %Modificamos input.in
    S = strrep(Finp,'Forward/Source999',['AdjointStn45/' DirSource]);              % MODIFICAR PATH
    S = strrep(S,'srfh','green');
    S = strrep(S,'number_of_point_sources = 1',['number_of_point_sources = ' num2str(1)]); % MODIFICAR PATH  number_of_point_sources = 1
    S = strrep(S,'number_output_planes        = 4',...
        ['number_output_planes        = ' num2str(size(Planes,1)+1)]);
    
    fileID = fopen('input.in','w');
    fprintf(fileID,S,'*write');
    fclose(fileID);
    clear S fileID;
    
    %Mesh of the domain
    for i=1:size(Planes,1)
        system(['cp ' PathAux Planes{i} ' .'])
    end
    
    %Creamos archivos de stationcoords
    fileID = fopen('planeinputcoords.0','w');
    fprintf(fileID,'%i\t',[0 1 1]);
    fprintf(fileID,'\n');
    fprintf(fileID,'%6.2f\t',[Eqk(iS).Loc.X Eqk(iS).Loc.Y Eqk(iS).Loc.Depth]);
    fclose(fileID);
    
    cd(TD);
end

