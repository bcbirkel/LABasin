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
PathSrc = '/home/dude/alanjuar/Documents/S20100707/Simulation/Forward/';

%% UserModification: Path to the Simulation directory.
%  Here we are going to create input directories and files
PathSim = '/home/dude/alanjuar/Documents/S20100707/Simulation/Adjoint/';

%% UserModification: Specify the desired minimum and maximum
% frequencies for waveform filtering. (it is a super good idea these match with
% the hercules simulation)
Mo= 4.0654e+20; % Seismic Moment from Hercules
FMax=1/6; %Upper Frequecy
FMin=1/50;  %Lower frequency
dt=0.25; 
t=0:dt:300;
iS=1;   %Source number
nE=16;   %Station number (depends on planeinput.0)
iP=2;   %Component (usually 1=x=lat, 2=y=lon, 3=z=depth)

%% UserModification: Path to the AuxFiles directory. Here we save or read
%  things that previous codes saved, for example the database.
PathAux = '/home/dude/alanjuar/Documents/S20100707/Codes/AuxFiles/';
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

SynF = bp_bu_co(Syn(iP).Signal,FMin,FMax,1/dt,4,2)';
ObsF = bp_bu_co(Obs(iP).Signal,FMin,FMax,1/dt,4,2)';

% Atomic Decomposition
[BASE,CTF] = atom_decm(SynF,dt,.5);
IsoF = sum(BASE(1:3,:),1)';

%%% Make Seismogram Decomposition %%%
subplot(3,2,1)
plot(t,ObsF,'k','linewidth',2), hold on
plot(t,SynF,'r','linewidth',2),
plot(t,IsoF,'-.','linewidth',1.5),
legend('Data U(t)','Synthetic S(t)','Isolation Filter I(t)','location','southeast')
xlim([0 t(end)]), grid minor, title('Data')
xlabel('Time [s]') , ylabel('Velocity [m/s]') 
xlim([50 250])

% Compute GSDF y correlations
Frec=0.025:0.025:0.15;

% Crosscorrelation and autocorrelation
no=length(SynF);
tt(:,1)=-dt*(no-1):dt:dt*(no-1);

ObsC(:,1)=xcorr(1e10*IsoF,1e10*ObsF,'biased');
SynC(:,1)=xcorr(1e10*IsoF,1e10*SynF,'biased');

% Pick the maximum of the correlation
Tw=2/Frec(ceil(length(Frec)/2));
sw=2*pi*0.72/Tw;    % Sigma w

% Find peaks of the crosscorrelation
[~,locs] = findpeaks(SynC);
[~,In]=min(abs(tt(locs)));
InS=locs(In);


[~,locs] = findpeaks(ObsC);
[~,In]=min(abs(tt(locs)-tt(InS)));
InO=locs(In);

% Windowing
WO=exp(-(0.5*sw^2)*(tt-tt(InO)).^2).*ObsC;
WS=exp(-(0.5*sw^2)*(tt-tt(InS)).^2).*SynC;  


subplot(3,2,2)
plot(tt,ObsC,'k','linewidth',2), hold on
plot(tt,SynC,'r','linewidth',2),
plot(tt,max(SynC)*exp(-(0.5*sw^2)*(tt-tt(InS)).^2),'-.','linewidth',2),
legend('C_I_U','C_I_S','Window W(t)','location','southeast')
grid minor, title('Crosscorrelograms'), xlim([-100 100])
xlabel('Lag Time [s]') , ylabel('Correlation Amplitude') 

subplot(3,2,3)
plot(tt,WO,'k','linewidth',2), hold on
plot(tt,WS,'r','linewidth',2),
plot(tt,max(SynC)*exp(-(0.5*sw^2)*(tt-tt(InS)).^2),'-.','linewidth',2),
legend('WC_I_U','WC_I_S','Window W','location','southeast')
grid minor, title('Windowed Crosscorrelograms'), xlim([-100 100])
xlabel('Lag Time [s]') , ylabel('Correlation Amplitude') 

% Parameters for "bootstraping"
si=0.1*2*pi*Frec(ceil(length(Frec)/2));   % Sigma i
dOR=computebandfftfilter_gauss(WO',dt,Frec(ceil(length(Frec)/2)),si);
dSR=computebandfftfilter_gauss(WS',dt,Frec(ceil(length(Frec)/2)),si);
[~,InOR]=max(dOR); % Maximo de las correlaciones
[~,InSR]=max(dSR);

% Compute for all frequencies
for nF=1:length(Frec)
    
    si=0.1*2*pi*Frec(nF);   % Sigma i
    % Crosscorrelagram and Autocorrelagram filtering
    dO=computebandfftfilter_gauss(WO',dt,Frec(nF),0.5*si);
    dS=computebandfftfilter_gauss(WS',dt,Frec(nF),0.5*si);

    % Let's do bootstraping
    [~,InO]=max(dO); 
    [~,InS]=max(dS);
    
    BS=1;
    while BS==1
        if (tt(InO) < tt(InOR)+0.5/Frec(nF)) && (tt(InO) >= tt(InOR)-0.5/Frec(nF))
            BS=0;
            break;
        elseif tt(InO) >= tt(InOR)+0.45/Frec(nF)
                InO=InO-round(1/Frec(nF)/dt);
        elseif tt(InO) < tt(InOR)-0.45/Frec(nF)
                InO=InO+round(1/Frec(nF)/dt);
        end
    end
    
    BS=1;
    while BS==1
        if (tt(InS) < tt(InSR)+0.5/Frec(nF)) && (tt(InS) >= tt(InSR)-0.5/Frec(nF))
            BS=0;
            break;
        elseif tt(InS) >= tt(InSR)+0.5/Frec(nF)
                InS=InS-round(1/Frec(nF)/dt);
        elseif tt(InS) < tt(InSR)-0.5/Frec(nF)
                InS=InS+round(1/Frec(nF)/dt);
        end
    end

    % Ajuste a la gaussiana (forma lenta)
    Eqn = 'a*exp(-0.5*b^2*(x-c)^2)*cos(d*(x-e))';

    GaS = fit(tt,dS,Eqn,...
        'Lower',[0.5*max(dS) 0.5*si tt(InS)-0.25/Frec(nF) 0.5*2*pi*Frec(nF) tt(InS)-0.25/Frec(nF)],...
        'Upper',[1.5*max(dS) 2.0*si tt(InS)+0.25/Frec(nF) 2.0*2*pi*Frec(nF) tt(InS)+0.25/Frec(nF)],...
        'StartPoint',[max(dS) si tt(InS) 2*pi*Frec(nF) tt(InS)]);

    GaO = fit(tt,dO,Eqn,...
        'Lower',[0.5*max(dO) 0.5*si tt(InO)-0.25/Frec(nF) 0.5*2*pi*Frec(nF) tt(InO)-0.25/Frec(nF)],...
        'Upper',[1.5*max(dO) 1.5*si tt(InO)+0.25/Frec(nF) 1.5*2*pi*Frec(nF) tt(InO)+0.25/Frec(nF)],...
        'StartPoint',[max(dO) si tt(InO) 2*pi*Frec(nF) tt(InO)]);

    
    %Parametros de GSDF
    w0=2*pi/(tt(end));               % Frecuencias
    wN=2*pi/(2*dt);
    w(:,1)=-wN:w0:wN;
    
    % Parametros:
    wi=2*pi*Frec(nF);
    
    wP=((si^2)*w+(sw^2)*wi)/(sw^2+si^2);
    wPP=((si^2)*w-(sw^2)*wi)/(sw^2+si^2);
    siP=((si^2)*(sw^2)/(sw^2+si^2))^0.5;
    
    tff=conj(fftshift(fft(IsoF)))*1/no;
    
%     % Check the GaS.a factor
    IW=(siP/(sw*GaS.a))*exp(-0.5*(w-2*pi*Frec(nF)).^2/(sw^2+si^2)).*tff./wP+...
        (siP/(sw*GaS.a))*exp(-0.5*(w+2*pi*Frec(nF)).^2/(sw^2+si^2)).*tff./wPP;
        
    IW(1:floor(size(IW)/2))=0*IW(1:floor(size(IW)/2));
    
    itff=ifft(fftshift(no*IW));
    
    % Save the GSDF measurements
    GSDF(nF,1)=Frec(nF);
    GSDF(nF,2)=GaO.e-GaS.e;	% delta_P
    GSDF(nF,3)=GaO.c-GaS.c;	% delta_G
    GSDF(nF,4)=-log(GaO.a/GaS.a)/GaS.d;  % Amplitud
    GSDF(nF,5)=-(GaO.d-GaS.d)/GaS.b^2;  % Amplitud
    
    %Save the seismogram perturbation kernel
    itff(isnan(itff))=0;
    Jpq(nF,:)=itff;
    
    % Save the Gaussian approximation parameters
    PO(nF,:)=[GaO.a GaO.b GaO.c GaO.d GaO.e];
    PS(nF,:)=[GaS.a GaS.b GaS.c GaS.d GaS.e];
    
    Fit = GaS.a*exp(-0.5*GaS.b^2*(tt-GaS.c).^2).*cos(GaS.d*(tt-GaS.e));
    if mod(nF,2) == 0
        subplot(3,2,4)
        plot(tt,0.5*nF+0.75*dS/max(dS),'k','linewidth',2), hold on
        plot(tt,0.5*nF+0.75*Fit/max(Fit),'-.r','linewidth',2), hold on
        text(-10,0.5*nF+0.1,[num2str(Frec(nF)) ' Hz)'])
        legend('F_iWC_I_U','F_iWC_I_S','location','southeast')
        xlabel('Lag Time [s]') ,
    end
    
    clear GaO GaS InO InS;
end

grid minor; xlim([-100 100])
title('Filtered Windowed Crosscorrelograms')

subplot(3,2,5)
plot(Frec,GSDF(:,2),'-*','linewidth',2), hold on
plot(Frec,GSDF(:,3),'-v','linewidth',2), hold on
plot(Frec,GSDF(:,4),'-s','linewidth',2), hold on
plot(Frec,GSDF(:,5),'-o','linewidth',2), hold on
xlabel('Frequency [Hz]'), grid minor, ylabel('\delta_t_g_q_a')
title('GSDF Measurements'), xlim([0 0.175])
legend('\delta\tau_t','\delta\tau_g','\delta\tau_q','\delta\tau_a')




