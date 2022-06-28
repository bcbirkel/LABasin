% Convert Hercules output to SAC files
% Written by Brianna Birkel 9/2021
clear all; close all; clc

%% UserModification: Path to the Source directory. In this directory are
% saved the outputs from the Hercules Simulation
PathSrc = '/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/runSim/Simulation/Forward/';

%% UserModification: Path to the Resources directory (catalogs, stations...)
PathRes='/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/runSim/Resources/';

%% UserModification: Path to the AuxFiles directory. Here we save or read
%  things that previous codes saved, for example the database.
PathAux = '/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/runSim/Codes/AuxFiles/';

% Loading Earthquakes and Stations for the Simulation
load([PathAux 'EarthquakeDB.mat']);

mkdir('/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/runSim/SACfiles')
%% UserModification: Specify the desired minimum and maximum
% frequencies for waveform filtering. (it is a super good idea these match with
% the hercules simulation)
FMax=1/4;
FMin=1/10;
Mo= 3.131630e+30; % Seismic Moment from Hercules MUST CHANGE BASED ON SLURM OUT!!
realMo = 6.29e+26;
Dt=0.05;
T=0:Dt:200;
iS=1;           % Number of Source to plot (SourceN)
% NumEst=[47    39    29    24    11     8    59    48    32    56    69    65 ...
%         12    51    20    43    25    60    30    28    17    64    18    22 ...
%         62    57    50    46    26    58    42    68    34    23    52    63 ...
%         31    36    45    37    66    33    55    40    44    67    38    49 ...
%         27    61    35    41    19     5    70    53    54     4    15    71];
NumEst=1:120;

%% No modification needed in this secction 16
% Loading information of the box simulation
FileID = fopen([PathAux 'XYtoLonLat.txt']);

Zone = fscanf(FileID, '%f', 1); Az = fscanf(FileID, '%f', 1);
X0_lat = fscanf(FileID, '%f', 1); Y0_lon = fscanf(FileID, '%f', 1);
XLim(1) = fscanf(FileID, '%f', 1); XLim(2) = fscanf(FileID, '%f', 1);
YLim(1) = fscanf(FileID, '%f', 1); YLim(2) = fscanf(FileID, '%f', 1);
ZLim(1) = fscanf(FileID, '%f', 1); ZLim(2) = fscanf(FileID, '%f', 1);
Corner(1) = fscanf(FileID, '%f', 1); Corner(2) = fscanf(FileID, '%f', 1);
fclose(FileID); An=abs(Az);

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


%% No modification needed in this secction
% Ploting recording at each station. Open plane
PathPlane=[PathSrc 'Source1/output/planes/planedisplacements.0'];
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
TimeSteps = round(T(end)/(Steps*DeltaT)+1);

% Reading Components for all stations
for nM=1:TimeSteps
    for iE=1:PointsStrike
        SynS(iE).X(nM) = fread(Plane, 1,'double');
        SynS(iE).Y(nM) = fread(Plane, 1,'double');
        SynS(iE).Z(nM) = fread(Plane, 1,'double');
    end
end
fclose(Plane);

% Rotation componets form LonLat to XY and plot
NF=1; % Number of figure
NS=1; %Subplot number

for nE=NumEst

    Factor = (Eqk(iS).Src.Mo*10^Eqk(iS).Src.Exp)/Mo;
%     Factor = realMo/Mo; 
    
    SX = Factor*(SynS(nE).X)
    SY = Factor*(SynS(nE).Y);
    SZ = -Factor*(SynS(nE).Z);

    %Rotate X,Y to R,T
    AnR=atan2(Eqk(iS).Stn(nE).X-Eqk(iS).Loc.X,Eqk(iS).Stn(nE).Y-Eqk(iS).Loc.Y)*180/pi;

    Syn(1).Signal = cosd(AnR)*SY+sind(AnR)*SX;
    Syn(2).Signal = -sind(AnR)*SY+cosd(AnR)*SX;
    Syn(3).Signal = SZ;
    
    clear SX SY SZ;
    
    t25 = tukeywin(length(Syn(1).Signal),0.05);
    CMP = ['-R';'-T';'-Z'];
     [Am,In]=max(max(abs([Syn(1).Signal; Syn(2).Signal; Syn(3).Signal])));
%     if Eqk(iS).Stn(nE).Num < 154
%         t0=max([40 min([5*floor(T(In)/5) 200-80])]);
%     end
%     if Eqk(iS).Stn(nE).Num > 154
%         t0=max([40 min([5*floor(T_iris(In)/5) 200-80])]);
%     end

    cd '/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/runSim/SACfiles'
    eventtime = datenum(Eqk(iS).Loc.Year,Eqk(iS).Loc.Month,Eqk(iS).Loc.Day,Eqk(iS).Loc.Hour,Eqk(iS).Loc.Min,Eqk(iS).Loc.Sec);
    mksac([Eqk(iS).Stn(nE).Name '_R.SAC'],Syn(1).Signal,eventtime,'DELTA',Dt,'KSTNM',[Eqk(iS).Stn(nE).Name],'KCMPNM','R');
    mksac([Eqk(iS).Stn(nE).Name '_T.SAC'],Syn(2).Signal,eventtime,'DELTA',Dt,'KSTNM',[Eqk(iS).Stn(nE).Name],'KCMPNM','T');
    mksac([Eqk(iS).Stn(nE).Name '_Z.SAC'],Syn(3).Signal,eventtime,'DELTA',Dt,'KSTNM',[Eqk(iS).Stn(nE).Name],'KCMPNM','Z');
    
% %% PLOTTING 
%     cd '/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/runSim/Codes/Figures'
%     figure(NF);
%     Error(nE)=0;
%     Eqk(iS).Stn(nE)
%     for iP=1:3
%         Sint=bp_bu_co(t25'.*Syn(iP).Signal,FMin,FMax,1/Dt,4,2);
%         subplot(3,3,NS+3*(iP-1))
%         plot(T,Sint,'r','linewidth',1.5),
%         xlabel 'Time [s]', 
%         xlim([t0-40 t0+80]), ylim(Am*[-1 1])
%         ylabel('Vel [m/s]')
%         text(t0+2.5-40,0.25*Am,[Eqk(iS).Stn(nE).Name  CMP(iP,:)])
%     end
%     NS=NS+1;
%     if NS>3
%         NS=1; NF=NF+1;
%         print(['Seismograms' num2str(NF-1) '.pdf'],'-dpdf','-bestfit')
%     end
%     clear Syn;
%     cd '/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/'
end
