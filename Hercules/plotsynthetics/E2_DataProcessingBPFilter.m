% Data processing using the S-Transform

% Written by Alan Juarez (UNAM-USC) Nov, 2018.
clear; close all; clc

%Paths to the outputs from the Hercules Simulations
PathSrc1 = '/home/scec-00/alanjuar/CoachellaValley/Simulation2/Forward/Source1/';
PathSrc2 = '/home/scec-00/alanjuar/CoachellaValley/Simulation2/Forward/Source2/';

%Specify the desired minimum and maximum frequencies for waveform filtering
FMax=1/4;
FMin=1/100;

% Load data from Simulation1
Plane1 = fopen([PathSrc1 '/output/planes/planedisplacements.0']);

Aux=  fread(Plane1, 8,'double');
BoxCorners=reshape(Aux,4,2);
LengthY=  fread(Plane1, 1,'double');
LengthX=  fread(Plane1, 1,'double');
PointsDip = fread(Plane1, 1,'int');
PointsStrike    = fread(Plane1, 1,'int');
DeltaT  = fread(Plane1, 1,'double');
Steps      = fread(Plane1, 1,'int');
TimeSteps = round(300/(Steps*DeltaT)+1);

for n = 1:TimeSteps
    for r = 1:PointsStrike
        SynS1(r).X(n) = fread(Plane1, 1,'double');
        SynS1(r).Y(n) = fread(Plane1, 1,'double');
        SynS1(r).Z(n) = fread(Plane1, 1,'double');
    end
end
fclose(Plane1);

% Load data from Simulation1
Plane2 = fopen([PathSrc2 '/output/planes/planedisplacements.0']);

Aux=  fread(Plane2, 8,'double');
BoxCorners=reshape(Aux,4,2);
LengthY=  fread(Plane2, 1,'double');
LengthX=  fread(Plane2, 1,'double');
PointsDip = fread(Plane2, 1,'int');
PointsStrike    = fread(Plane2, 1,'int');
DeltaT  = fread(Plane2, 1,'double');
Steps      = fread(Plane2, 1,'int');
TimeSteps = round(300/(Steps*DeltaT)+1);

for n=1:TimeSteps
    for r=1:PointsStrike
        SynS2(r).X(n) = fread(Plane2, 1,'double');
        SynS2(r).Y(n) = fread(Plane2, 1,'double');
        SynS2(r).Z(n) = fread(Plane2, 1,'double');
    end
end
fclose(Plane2);

% Simulation parameters
Dt = (Steps*DeltaT);
T = 0:Dt:300;
nt = length(T);

% frequencies needed
fN = 1/(2*Dt);
f0 = 1/(T(end)-T(1));
Fr = 0:f0:fN;

Mo = 3.553074e+30; % Seismic Moment from Hercules
load([ './AuxFiles/EqkStnDB.mat']);
Factor = (Eqk(1).Src.Mo*10^Eqk(1).Src.Exp)/Mo;

% Select some stations
iS = [88 111 106 29 34 100 57 68 37 11 74]; %[34 37 111 11 100 106 74 57];%
  
% Create Periods with 4*Lambda < Dist
MinP = 5;
MaxP = 50;

np=25;
for nf=1:np
    Periods(nf)=MinP+exp((log(MaxP-MinP+1))/(np-1)*(nf-1))-1;
end
Freqs = 1./Periods;


% Select appropriate velocities
MinC=0.1e3;
MaxC=1.5e3;
for nf=1:np
    C(nf)=MinC+exp((log(MaxC-MinC+1))/(np-1)*(nf-1))-1;
end


  cd FiguresSim2Auto
  
for r = iS
    % Correct Hercules Moment
    SX1 = Factor*(SynS1(r).X);
    SY1 = Factor*(SynS1(r).Y);
    SZ1 = -Factor*(SynS1(r).Z);
    
    SX2 = Factor*(SynS2(r).X);
    SY2 = Factor*(SynS2(r).Y);
    SZ2 = -Factor*(SynS2(r).Z);
    
    % Rotate X,Y to R,T
    AnR = atan2(Stn(r).X-Eqk.Loc.X,Stn(r).Y-Eqk.Loc.Y)*180/pi;
    
    Syn1(1,:) = cosd(AnR)*SY1+sind(AnR)*SX1;
    Syn1(2,:) = -sind(AnR)*SY1+cosd(AnR)*SX1;
    Syn1(3,:) = SZ1;
    
    Syn2(1,:) = cosd(AnR)*SY2+sind(AnR)*SX2;
    Syn2(2,:) = -sind(AnR)*SY2+cosd(AnR)*SX2;
    Syn2(3,:) = SZ2;
    
    % Extra parameters
    t25 = tukeywin(length(Syn1(1,:)),0.25);
    CMP = ['R';'T';'Z'];
    [Am,~] = max(abs(Syn1(:)));
    
    % Distance
    Dist = 0.001*norm([Eqk.Loc.X Eqk.Loc.Y]-[Stn(r).X Stn(r).Y]);
    
    figure;
    for c = 1:3
        % Filtering
        S1 = bp_bu_co(t25'.*Syn1(c,:),FMin,FMax,1/Dt,4,2);
        S2 = bp_bu_co(t25'.*Syn2(c,:),FMin,FMax,1/Dt,4,2);
        
        %fft and frequency axis
        Sfft1 = fft(S1);
        Sfft2 = fft(S2);
        H1 = Sfft1(1:length(Fr));
        H2 = Sfft2(1:length(Fr));

        % Generate narrow band gaussian filtered data
        for nf = 1:length(Freqs)
            Gausf = exp(-C(nf)*(Fr-Freqs(nf)).^2./2./Freqs(nf));
            F1 = H1.*Gausf;
            F1(end+1:nt) = 0;
            nband1(nf,:) = ifft(F1);
            A1(nf) = find(abs(nband1(nf,:))==max(abs(nband1(nf,:))));

            F2 = H2.*Gausf;
            F2(end+1:nt) = 0;
            nband2(nf,:) = ifft(F2);
            A2(nf) = find(abs(nband2(nf,:))==max(abs(nband2(nf,:))));
        end
        
        C1(c,:) = smooth(Dist./T(A1));
        C2(c,:) = smooth(Dist./T(A2));
        dC(c,:) = 1e2*(C2(c,:)-C1(c,:))./C1(c,:);
        
        subplot(3,3,c)
        plot(T,S1,'color',[0 0 0],'linewidth',1.5), hold on
        plot(T,S2,'color',[0.635 0.078 0.184],'linewidth',1),
        xlabel 'Time (s)', ylabel([CMP(c,:) ' (m)']) %legend('Observado','Sintetico'), 
        xlim([0 200]), ylim(Am*[-1 1]), axis square
        title([Stn(r).Name])
        
        subplot(3,3,c+3)
        plot(Periods,C1(c,:),'.','color',[0 0 0],'linewidth',1.5); hold on
        plot(Periods,C2(c,:),'.','color',[0.635 0.078 0.184],'linewidth',1);
        xlabel 'Period (s)', ylabel 'c (km/s)';
        ylim([0.2 4.2]), xlim([Periods(1) Periods(end)]), axis square
        
%         subplot(4,3,c+6)
%         semilogy(Periods(dC(c,:) > 0), dC(c,dC(c,:) > 0),'.','color',[0  0.447  0.741],'linewidth',2); hold on
%         semilogy(Periods(dC(c,:) <= 0), -dC(c,dC(c,:) <= 0),'.','color',[0.635 0.078 0.184],'linewidth',2); hold on
%         xlabel 'Period (s)', ylabel '|\Deltac/c| (%)'; 
%         ylim([1e-2 1e2]),xlim([Periods(1) Periods(end)]), axis square
        
        subplot(3,3,c+6)
        plot(Periods, dC(c,:),'.','color',[0.635 0.078 0.184],'linewidth',2); hold on
        xlabel 'Period (s)', ylabel '(\Deltac)/c (%)'; 
        ylim(10*[-1 1]), xlim([Periods(1) Periods(end)]), axis square
        
    end
   
    save (['Sim2DispCurves' num2str(r) '.mat'],'Periods','C1','C2','dC');
    print(['Sim2DispCurves' num2str(r) '.pdf'],'-dpdf','-bestfit');
    
    print(['Sim2DispCurves' num2str(r) '.png'],'-dpng')
    close;
    clear C1 C2 dC;
end
cd ..
