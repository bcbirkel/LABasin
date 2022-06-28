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
PathSrc = '/home/scec-00/alanjuar/S20100707/Simulation/Forward/';

%% UserModification: Path to the Simulation directory.
%  Here we are going to create input directories and files
PathSim = '/home/scec-00/alanjuar/S20100707/Simulation/Adjoint/';

%% UserModification: Specify the desired minimum and maximum
% frequencies for waveform filtering. (it is a super good idea these match with
% the hercules simulation)
FMax=1/4;
FMin=1/30;
Mo= 4.065432e+30; % Seismic Moment from Hercules
dt=0.25;
t=0:dt:200;
iS=1;   %Source number
nE=18;  %Station number (depends on planeinput.0)
iP=2;   %Component (usually 1=x=lat, 2=y=lon, 3=z=depth)

%% UserModification: Path to the AuxFiles directory. Here we save or read
%  things that previous codes saved, for example the database.
PathAux = '/home/scec-00/alanjuar/S20100707/Codes/AuxFiles/';
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
SZ = -Factor*SynS(nE).Z;

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

%% Atomic Decomposition
figure;
[BASE,CTF] = atom_decm(SynF,dt,FMin,FMax,1/5,.5);


%%% Make Seismogram Decomposition %%%
subplot(2,2,1)
plot(t,SynF,'k','linewidth',2), hold on

subplot(2,2,1)
plot(t,sum(BASE,1),'r','linewidth',1.5),
legend('Data','\SigmaPhases','location','southeast')
xlim([0 t(end)]), ylim(max(abs(SynF))*1.5*[-1 1]),
title('Reconstruction')
xlabel('Time [s]') , ylabel('Displacement [m]') 

MMax=10;%size(BASE,1);

subplot(2,2,[2 4])
for m = 1:MMax 
    plot(t,m+1.5*BASE(m,:)/max(max(abs(BASE))),'r','linewidth',1.5), hold on
end
xlim([0 t(end)]), ylim([0 m+1]), 
xlabel ('Time [s]'), ylabel ('Phases')
title('Dictionary'), 

%% Compute GSDF y correlations
Frec=[1/20 1/15 1/10 1/7.5 1/5];
for nF=1:length(Frec)
    [FWC.GSDF(nF,:),FWC.Jpq(nF,:),FWC.PS(nF,:),~]=gsdf(1e10*SynF,1e10*SynF,1e10*SynF,dt,Frec(nF));
end
% Uncomment if you use a factor in the GSDF function to avoid small numbers
FWC.Jpq=1e10*FWC.Jpq;
FWC.PS(:,1)=FWC.PS(:,1)/1e20;

for m=1:MMax
    for n=1:MMax
        for nF=1:length(Frec)
            [FWCmn(m,n).GSDF(nF,:),FWCmn(m,n).Jpq(nF,:),FWCmn(m,n).PS(nF,:),~] = ...
                gsdf(1e10*BASE(n,:)',1e10*SynF,1e10*BASE(m,:)',dt,Frec(nF));
        end
        % Uncomment if you use a factor in the GSDF function to avoid small numbers
        FWCmn(m,n).Jpq=1e10*FWCmn(m,n).Jpq;
        FWCmn(m,n).PS(:,1)=FWCmn(m,n).PS(:,1)/1e20;
    end
end

%Slip function
FWC.ASF=sum(imag(FWC.Jpq),1);

mkdir(PathSim)
for m=1:MMax
    for nF=1:length(Frec)
        %Path to simulation and forlders%
        cd (PathSim)
        
        DirSource= ['Phase' num2str(m) '_Freq' num2str(nF)];
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
        fprintf(fileID,'%i\t',length(FWCmn(m,m).Jpq(nF,:)));
        fprintf(fileID,'%i\t',0);
        fprintf(fileID,'%6.4f\t',dt);
        fprintf(fileID,'%e\t',imag(FWCmn(m,m).Jpq(nF,:)));
        fprintf(fileID,'\n');
        fclose(fileID);
        clear fileID;
        
        %Modificamos input.in
        S = strrep(Finp,'Forward/Source999',['Adjoint/' DirSource]);              % MODIFICAR PATH
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
end

%% Correlation Reconstruction
tt(:,1)=-dt*(length(t)-1):dt:dt*(length(t)-1);

for nF=1:length(Frec)
    FWC.FWC(nF,:)=FWC.PS(nF,1)*exp(-0.5*FWC.PS(nF,2)^2*(tt-FWC.PS(nF,3)).^2)...
        .*cos(FWC.PS(nF,4)*(tt-FWC.PS(nF,5)));
    for m=1:MMax%size(BASE,1)
        for n=1:MMax%size(BASE,1)
            FWCmn(m,n).FWC(nF,:)=FWCmn(m,n).PS(nF,1)*exp(-0.5*FWCmn(m,n).PS(nF,2)^2*(tt-FWCmn(m,n).PS(nF,3)).^2)...
                .*cos(FWCmn(m,n).PS(nF,4)*(tt-FWCmn(m,n).PS(nF,5)));
        end
    end
end

%% Interference factors
for nF=1:length(Frec)
    [wf(nF,:),sf(nF,:)]=compute_wf_sf(SynF,SynF,dt,Frec(nF));
end

FWC.tp=FWC.PS(:,5);
FWC.tg=FWC.PS(:,3);
FWC.tq=-log(FWC.PS(:,1))./wf(:);
FWC.ta=-(FWC.PS(:,4)-wf(:))./(sf(:).^2);

for m=1:MMax
    for n=1:MMax
        for nF=1:length(Frec)
            FWCmn(m,n).tp(nF)=FWCmn(m,n).PS(nF,5);
            FWCmn(m,n).tg(nF)=FWCmn(m,n).PS(nF,3);
            FWCmn(m,n).tq(nF)=-log(FWCmn(m,n).PS(nF,1))/wf(nF);
            FWCmn(m,n).ta(nF)=-(FWCmn(m,n).PS(nF,4)-wf(nF))/sf(nF)^2;
            B(m,n,nF)=exp(-wf(nF)*(FWCmn(m,n).tq(nF)-FWC.tq(nF)))...
                *exp(-0.5*sf(nF)^2*(FWCmn(m,n).tg(nF)-FWC.tg(nF))^2);
            P(m,n,nF)=(wf(nF)-sf(nF)^2*FWCmn(m,n).ta(nF))*...
                (FWCmn(m,n).tp(nF)-FWC.tg(nF))-wf(nF)*(FWC.tp(nF)-FWC.tg(nF));
        end
    end
end


C=B.*cos(P);
CM = loadcmap('BlueWhiteOrangeRed.c3g');

for nF = 1:3:length(Frec)
    figure;
    for ia = 1:3
        if nF+ia-1 <= length(Frec)
            sp(nF+ia-1) = subplot(3,3,3*ia-2);
            imagesc(B(:,:,nF+ia-1).^0.5)
            colormap(sp(nF+ia-1),CM),
            xlabel('n'), ylabel('m'), axis equal
            xlim([0.5 length(B(:,:,nF+ia-1))+0.5]),
            c1(nF+ia-1) = colorbar('Limits',0.5*max(max(abs(B(:,:,nF+ia-1))))*[0 1]);
            title(['B_m_n(' num2str(Frec(nF+ia-1)) 'Hz)'],'fontsize',9),
            caxis(sp(nF+ia-1),0.5*max(max(abs(B(:,:,nF+ia-1))))*[-1 1])
            
            sp(nF+ia-1) = subplot(3,3,3*ia-1)
            imagesc(P(:,:,nF+ia-1))
            colormap(sp(nF+ia-1),CM),
            xlabel('n'), ylabel('m'), axis equal
            xlim([0.5 length(P(:,:,nF+ia-1))+0.5])
            c1(nF+ia-1) = colorbar('Limits',0.5*max(max(abs(P(:,:,nF+ia-1))))*[-1 1]);
            title(['P_m_n(' num2str(Frec(nF+ia-1)) 'Hz)'],'fontsize',9),
            caxis(sp(nF+ia-1),0.5*max(max(abs(P(:,:,nF+ia-1))))*[-1 1])
            
            sp(nF+ia-1) = subplot(3,3,3*ia-0)
            imagesc(abs(C(:,:,nF+ia-1)).^0.5)
            colormap(sp(nF+ia-1),CM),
            xlabel('n'), ylabel('m'), axis equal
            xlim([0.5 length(C(:,:,nF+ia-1))+0.5])
            c1(nF+ia-1) = colorbar('Limits',0.5*max(max(abs(C(:,:,nF+ia-1))))*[0 1]);
            title(['C_m_n(' num2str(Frec(nF+ia-1)) 'Hz)'],'fontsize',9),
            caxis(sp(nF+ia-1),0.5*max(max(abs(C(:,:,nF+ia-1))))*[-1 1])
        end
    end
end


ASF_SS=0*t;
for m=1:MMax
    for nF=1:length(Frec)
        ASF_SS=ASF_SS+C(m,m,nF)*imag(FWCmn(m,m).Jpq(nF,:));
    end
end

figure;
subplot(2,2,1)
plot(t,FWC.ASF,'k','linewidth',2), hold on
plot(t,ASF_SS,'r','linewidth',2), hold on
xlabel('Time [s]'), grid minor, legend('Syn','\SigmaC_m_nU_n')
grid, title('Seismogram Perturbation Kernel'), xlim([t(1) t(end)])

save('./AuxFiles/InterFactFrec.mat','FWC')
save('./AuxFiles/InterFactFrec.mat','FWCmn','-append')
save('./AuxFiles/InterFactFrec.mat','C','-append')
save('./AuxFiles/InterFactFrec.mat','BASE','-append')