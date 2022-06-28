% Full Waveform Inversion - Code No.4.
% A3_EarthquakeDatabase.m builds up the source-station-waveform database.
% This code looks at all events in the catalog created in previous codes
% and looks for the waveforms to put all together. The final database
% contains velocities.
% User needs to go througth code and set all the parameters specified with
% "UserModification" legend on the comments.
% Special care must be taken in the folder structure where catalogs,
% resources and databases are located or where they are going to be saved.


clear all; close all; clc

%% UserModification: Specify the desired time-sampling and minimum and maximum
% frequencies for waveform filtering. (it is a super good idea these match with
% the hercules simulation, specially the final sampling)
FinalSPS=4;
Duration=200;
FMax=1;
FMin=1/100;

%% UserModification: Path to the AuxFiles directory. Here we save or read
%  things that previous codes saved, for example the database.
PathAux = '/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/Codes/AuxFiles/';

%% UserModification: Path to the Resources directory (catalogs, stations...)
% Here we need acces to stations since in this folder we want to find the
% instrument response files.
PathRes='/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/ObsData/';

%% UserModification: Path to the waveform directory (here we must find two
% folders, one with the sac velocity files and other with the acceleration
% waveforms. Reading seisan files for acceleration is complicated because
% you need a station list. The function I use here works for the BMSM.
PathWF='/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/ObsData/';

%% No modification needed in this secction
% Starting to sort database and waveforms
load ([PathAux 'EqkStnDB.mat']);

Stn=Data.Stn';
Eqk=Data.Eqk';
clear Data;

%% No modification needed in this secction
% Loading information of the box simulation
% Rotation data. We are going to set the stations coordinates in the XY
% system.
FileID = fopen([PathAux 'XYtoLonLat.txt']);

Zone = fscanf(FileID, '%f', 1); Az = fscanf(FileID, '%f', 1);
X0_lat = fscanf(FileID, '%f', 1); Y0_lon = fscanf(FileID, '%f', 1);
XLim(1) = fscanf(FileID, '%f', 1); XLim(2) = fscanf(FileID, '%f', 1);
YLim(1) = fscanf(FileID, '%f', 1); YLim(2) = fscanf(FileID, '%f', 1);
ZLim(1) = fscanf(FileID, '%f', 1); ZLim(2) = fscanf(FileID, '%f', 1);
Corner(1) = fscanf(FileID, '%f', 1); Corner(2) = fscanf(FileID, '%f', 1);
fclose(FileID);

%% UserModification: Read Carefully
% In this secction we are going to look for all waveform files. The velocity 
% files should be in SAC format. 
% This code will look for the three components (N,E,Z) in the 
% same directory and it expect them to have generic names in the IRIS 
% format. The code also needs to read individual files for each component 
% for each station for each earthquake.
AVel=get_list_files([PathWF '/Velocity/'],'*.SAC');

%% No User Modification Needed in this Section

for iS=1:size(Eqk,1)
    % Look for the file to read. Possibly you have to edit format.
    ReadFile='2014.088.040942';
    
    % Look for velocity files, read them and do the same data processing to
    % make everything uniform
    Vel=strfind(AVel,ReadFile);
    iE=1; iC=[]; ARead=0;
    for iAv=1:size(Vel,1)
        if isempty(Vel{iAv})==0
            
            SF=rsac('little-endian',[PathWF '/Velocity/' char(AVel{iAv})]);
            u1.dep=SF(:,2); u1.head=SF(:,3); clear SF;
            ARead=ARead+1;
            ArS=char(AVel{iAv});
            IP=find(ArS=='.');
            File(ARead).Name=ArS(IP(1)+1:IP(2)-1);
            if ARead>1
                if strcmp(File(ARead).Name,File(ARead-1).Name)==0
                    iE=iE+1;
                end
            end
            
                if ArS(IP(4)-1)=='N'
                    iC=1;
                elseif ArS(IP(4)-1)=='E'
                    iC=2;
                elseif ArS(IP(4)-1)=='Z'
                    iC=3;
                end
                Eqk(iS).Stn(iE).Cmpt(iC).Name=[ArS(IP(1)+1:IP(2)-1) ArS(IP(3)+1:IP(4)-1)];
                Eqk(iS).Stn(iE).Name=ArS(IP(1)+1:IP(2)-1);
            
            if abs(u1.head(32)) >= 90
                for ii=1:size(Stn,1)
                    if strcmp(Eqk(iS).Stn(iE).Name,Stn(ii).Name)
                        Eqk(iS).Stn(iE).Lon=Stn(ii).Lon;
                        Eqk(iS).Stn(iE).Lat=Stn(ii).Lat;
                        Eqk(iS).Stn(iE).Elev=Stn(ii).Elev;
                        Eqk(iS).Stn(iE).Resp=Stn(ii).Resp;
                        Eqk(iS).Stn(iE).Num=ii;
                        break;
                    end
                end
            else
                Eqk(iS).Stn(iE).Lon=u1.head(33);
                Eqk(iS).Stn(iE).Lat=u1.head(32);
                Eqk(iS).Stn(iE).Elev=u1.head(34);
                for ii=1:size(Stn,1)
                    if strcmp(Eqk(iS).Stn(iE).Name,Stn(ii).Name)
                        Eqk(iS).Stn(iE).Num=ii;
                        break;
                    end
                end
            end
            Eqk(iS).Stn(iE).Cmpt(iC).SPS=FinalSPS;
            
            Sps=round(1/u1.head(1));
            
            Sint= rmean(rtrend(u1.dep))';
            u1.dep=Sint;
            clear Sint;
                    
            Avg=4*mean(abs(u1.dep(1:2*Sps)));
            if length(u1.dep)<Duration*Sps+1
                Signal=[u1.dep; Avg*rand([Duration*Sps+1-length(u1.dep) 1])-0.5*Avg];
            else
                Signal=u1.dep(1:Duration*Sps+1);
            end
            
            TW=tukeywin(length(Signal),0.01);
            Signal=bp_bu_co(TW.*Signal,FMin,FMax,Sps,2,2);
            clear TW a;
            
            if Sps <= FinalSPS
                Signal=interp(Signal,round(FinalSPS/Sps));
            else
                Signal=decimate(Signal,round(Sps/FinalSPS));
            end
            
            TW=tukeywin(length(Signal),0.01);
            if isempty(ArS(IP(2)+1:IP(3)-1))==1
                StnPos='--';
            else
                StnPos=ArS(IP(2)+1:IP(3)-1);
            end
                
            Signal=transfer(Signal,length(Signal),FinalSPS,...
                [PathRes 'Velocity/' 'SACPZ.' ArS(1:IP(2))  ...
                StnPos ArS(IP(3):IP(4)-1)],FMin,FMax);
            clear TW;
            
            Signal=cast(Signal,'double');
            
% %             %Uncoment this section id you preffer displacements
% %             Signal=cumtrapz(Signal)/FinalSPS;
% %             Signal = Signal-conv(sgolayfilt(Signal,2,a),(1/a)*ones([1 a]),'same');
            
            Avg=4*mean(abs(rtrend(rmean(Signal(1:5*FinalSPS)))));
            if length(Signal)<Duration*FinalSPS+1
                Signal=[Signal; Avg*rand([Duration*FinalSPS+1-length(Signal) 1])-0.5*Avg];
            else
                Signal(Duration*FinalSPS+2:end)=[];
            end
            
            TW=tukeywin(length(Signal),0.01)';
            Signal=bp_bu_co(TW.*Signal,FMin,FMax,FinalSPS,4,2);
            
            u1.dep=[];
            u1.dep=Signal;
            Eqk(iS).Stn(iE).Cmpt(iC).Signal=u1.dep;
            clear Signal Sps Avg u1 TW;
            
            
        end
    end
    
    % Select traces with signal to noise ratio greater than 8. You can
    % change this of course.
    CleanS=[]; nS=1; CMP=['BHN';'BHE';'BHZ';];
    for iA=1:size(Eqk(iS).Stn,2)
        if isempty (Eqk(iS).Stn(iA).Num)==1
            CleanS(nS)=iA;
            nS=nS+1;
        else
            for iC=1:size(Eqk(iS).Stn(iA).Cmpt,2)
                if isempty(Eqk(iS).Stn(iA).Cmpt(iC).Name)==1
                    Eqk(iS).Stn(iA).Cmpt(iC).Name=[Eqk(iS).Stn(iA).Name CMP(iC,:)];
                    Eqk(iS).Stn(iA).Cmpt(iC).SPS=FinalSPS;
                    Eqk(iS).Stn(iA).Cmpt(iC).Signal=zeros([1 Duration*FinalSPS+1]);
                else
                    Signal=Eqk(iS).Stn(iA).Cmpt(iC).Signal;
                    Dt=1/Eqk(iS).Stn(iA).Cmpt(iC).SPS;
                    
                    Dist=vdist(Eqk(iS).Loc.Lat,Eqk(iS).Loc.Lon,...
                        Eqk(iS).Stn(iA).Lat,Eqk(iS).Stn(iA).Lon)/1000;
                    tauPP=tauptime('mod','iasp91','dep',...
                        Eqk(iS).Loc.Depth*0.001,'km',Dist);
                    
                    uP=round(tauPP(1).time*Eqk(iS).Stn(iA).Cmpt(iC).SPS);
                    
                    if uP>length(Signal)
                        uP=length(Signal);
                    end
                    [Ener1,~,~,~]=movint(Signal(1:uP),Dt,.01,.99);
                    if 4*uP > length(Signal)
                        [Ener2,~,~,~]=movint(Signal(uP:end),Dt,.01,.99);
                    else
                        [Ener2,~,~,~]=movint(Signal(uP:2*uP),Dt,.01,.99);
                    end
                    
                    if (3*Ener1 < Ener2) %&& norm(Signal) < 1
%                         figure(iA)
%                         subplot(3,1,iC)
%                         Tiempo=0:Dt:Dt*(length(Signal)-1);
%                         plot(Tiempo,Signal), hold on
                    else
                        Eqk(iS).Stn(iA).Cmpt(iC).Name=[Eqk(iS).Stn(iA).Name CMP(iC,:)];
                        Eqk(iS).Stn(iA).Cmpt(iC).SPS=FinalSPS;
                        Eqk(iS).Stn(iA).Cmpt(iC).Signal=zeros([1 FinalSPS*Duration+1]);
                    end
                    clear Signal Dist uP Ener1 Ener2 tauPP Dt;
                end
            end
        end
    end
    if isempty(CleanS)==0
        Eqk(iS).Stn(CleanS)=[];
    end
        clear CleanS nS;
    
    % Set stations and sources in the XY coordinate system. Rotation
    % required, Maybe I should rotate to Radial Transverse but never mind.
    for iA=1:size(Eqk(iS).Stn,2)
        [Y,X] = ll2utm(Eqk(iS).Stn(iA).Lat,Eqk(iS).Stn(iA).Lon,Zone);
        CoorXY=[[cosd(Az) -sind(Az); sind(Az) cosd(Az)]*[Y-Y0_lon; X-X0_lat]]';
        
        Eqk(iS).Stn(iA).X=CoorXY(2);
        Eqk(iS).Stn(iA).Y=CoorXY(1);
        clear X Y CoorXY
    end
    [Y,X] = ll2utm(Eqk(iS).Loc.Lat,Eqk(iS).Loc.Lon,Zone);
    CoorXY=[[cosd(Az) -sind(Az); sind(Az) cosd(Az)]*[Y-Y0_lon; X-X0_lat]]';
    Eqk(iS).Loc.X=CoorXY(2);
    Eqk(iS).Loc.Y=CoorXY(1);
end

nB=1; Clean=[];
for iS=1:size(Eqk,1)
    if isempty(Eqk(iS).Stn) || size(Eqk(iS).Stn,2)==0
        Clean(nB)=iS;
        nB=nB+1;
    else
        CleanS=[]; nS=1;
        for iA=1:size(Eqk(iS).Stn,2)
            if size(Eqk(iS).Stn(iA),2)>2
                Eqk(iS).Stn(iA)
                if norm(Eqk(iS).Stn(iA).Cmpt(1).Signal)==0 ...
                        && norm(Eqk(iS).Stn(iA).Cmpt(2).Signal)==0 ...
                        && norm(Eqk(iS).Stn(iA).Cmpt(3).Signal)==0
                    CleanS(nS)=iA;
                    nS=nS+1;
                end
            end
        end
        if isempty(CleanS)==0
            Eqk(iS).Stn(CleanS)=[];
        end
        clear CleanS nS;
    end
end
Eqk(Clean)=[];

iF=1; Clean=[];
for iCS=1:size(Eqk,1)
    if size(Eqk(iCS).Stn,2)==0 || size(Eqk(iCS).Stn,2) < 6
        Clean(iF)=iCS;
        iF=iF+1;
    end
end
Eqk(Clean)=[];

% Finaly save the database with Earthquakes, Stations and Waveforms
clear Clean nB;
save([PathAux 'EarthquakeDB.mat'],'Eqk');

