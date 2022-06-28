% Data processing using the S-Transform

% Written by Alan Juarez (UNAM-USC) Nov, 2018.
clear; close all; clc

%Paths to the outputs from the Hercules Simulations
PathSrc1 = '/home/scec-00/alanjuar/CoachellaValley/Simulation2/Forward/Source1/';
PathSrc2 = '/home/scec-00/alanjuar/CoachellaValley/Simulation2/Forward/Source2/';

%Specify the desired minimum and maximum frequencies for waveform filtering
FMax=1/5;
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
Mo = 3.553074e+30; % Seismic Moment from Hercules
load([ './AuxFiles/EqkStnDB.mat']);
Factor = (Eqk(1).Src.Mo*10^Eqk(1).Src.Exp)/Mo;

% Select some stations
iS = 81;

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
    
    t25 = tukeywin(length(Syn1(1,:)),0.25);
    CMP = ['-R';'-T';'-Z'];
    [Am,~] = max(abs(Syn1(:)));
    
    figure;
    for c = 1:3
        S1 = bp_bu_co(t25'.*Syn1(c,:),FMin,FMax,1/Dt,4,2);
        S2 = bp_bu_co(t25'.*Syn2(c,:),FMin,FMax,1/Dt,4,2);
        
        subplot(3,1,c)
        plot(T,S1,'color',[0 0 0],'linewidth',2), hold on
        plot(T,S2,'color',[0.635 0.078 0.184],'linewidth',1.5),
        xlabel 'Time (s)', ylabel('Disp. (m)') %legend('Observado','Sintetico'), 
        xlim([0 300]), ylim(Am*[-1 1])
        text(T(1)+2.5,0.8*Am,[Stn(r).Name  CMP(c,:)])
    end
    
    % some data
    ns = length(T);
    ds1 = reshape(S1,1,ns);
    ds2 = reshape(S2,1,ns);

    % frequencies needed
    fN = 1/(2*Dt);
    f0 = 1/(T(end)-T(1));
    fr = 0:f0:fN;

    % S-transform of the original signal
    ST_s1 = stran(ds1,1);
    ST_s2 = stran(ds2,1);
    
    Dist = 0.001*norm([Eqk.Loc.X Eqk.Loc.Y]-[Stn(r).X Stn(r).Y]);
    IT1 = round(Dist/3/Dt);
    IT2 = round(Dist/1/Dt);
    
    for f = 1:length(fr)
        Freq(f) = fr(f);
        if fr(f) >= FMin && fr(f) <= FMax
            W = repmat(normpdf(fr,fr(f),f0)',1,length(T));

            IF = istransform(W.*ST_s1);
            [~,In] = max(IF(IT1:IT2));
            Tao1(f) = T(In+IT1);


            IF = istransform(W.*ST_s2);
            [~,In] = max(IF(IT1:IT2));
            Tao2(f) = T(In+IT1);
        else
            Tao1(f) = 0;
            Tao2(f) = 0;
        end
    end
    
    figure; 
    subplot(3,1,1)
    plot(T,ds1,'color',[0 0 0],'linewidth',2), hold on
    plot(T,ds2,'color',[0.635 0.078 0.184],'linewidth',1.5),
    xlabel 'Time (s)', ylabel('Disp. (m)'), legend('CVM-H','dVs'), 
    xlim([T(1) T(end)]), ylim(Am*[-1 1])
    text(T(1)+2.5,0.8*Am,[Stn(r).Name  CMP(c,:)])
        
    subplot(3,1,2) 
    surf(T,fr,abs(ST_s1).^2)
    shading interp, view (2), hold on , set(gca,'yscale','log'), 
    title ('Synthetic'), ylim([FMin FMax]), xlim([T(1) T(end)])
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    CM = loadcmap('BlueWhiteOrangeRed.c3g'); colormap(CM), 
    plot3(Tao1,Freq,1e10+0*Tao1,'.','color',[0 0 0],'markersize',8), hold on
    
    subplot(3,1,3) 
    surf(T,fr,abs(ST_s2).^2)
    shading interp, view (2), hold on , set(gca,'yscale','log'), 
    title ('Synthetic'), ylim([FMin FMax]), xlim([T(1) T(end)])
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    CM = loadcmap('BlueWhiteOrangeRed.c3g'); colormap(CM), 
    plot3(Tao2,Freq,1e10+0*Tao2,'.','color',[0.635 0.078 0.184],'markersize',8)
    
    figure;
    plot(1./Freq,Dist./Tao1,'.','color',[0 0 0],'markersize',10), hold on
    plot(1./Freq,Dist./Tao2,'.','color',[0.635 0.078 0.184],'markersize',10)
    xlabel('Period'), ylabel('C (m/s)')
    
end

