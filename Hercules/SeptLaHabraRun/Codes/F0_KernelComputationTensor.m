% Full Waveform Inversion - Code No F1.
% F1_KernelComputation.m computes the Full Waveform Inversion Kernel in all
% the domain reading the outputs from the Forward and Adjoint Hercules
% SImulations.
% Written by Alan Juarez (UNAM-USC) Dec, 2016.
clear; close all; clc;

%% UserModification: Path to the Source directory. In this directory are
% saved the Forward outputs from the Hercules Simulation
PathFor = '/home/scec-00/alanjuar/S20100707/Simulation/Forward/';

%% UserModification: Path to the Source directory. In this directory are
% saved the Adjoint outputs from the Hercules Simulation
PathAdj = '/home/scec-00/alanjuar/S20100707/Simulation/Adjoint/';

%% UserModification: Path to the AuxFiles directory. Here we save or read
%  things that previous codes saved, for example the database.
PathAux = '/home/scec-00/alanjuar/S20100707/Codes/AuxFiles/';

%% No User Modification Needed: Creating Adjoint Simulation Files
TD = pwd;
t = 200;  %Simulation Duration

% Estimate moment factor correction 
load([PathAux 'EarthquakeDB.mat']);
Mo= 4.065432e+30; % Seismic Moment from Hercules
if Mo == 0
    Factor = 1;
else
    Factor = (Eqk(1).Src.Mo*10^Eqk(1).Src.Exp)/Mo;
end
clear Eqk;
    
%Reading generic planeinputcoords file
Planes = get_list_files(PathAux, 'planeinputcoords.*');
for ii=1:3:size(Planes,1)
    PlaneInput = load([PathAux Planes{ii}]);
    NumP((ii+2)/3) = PlaneInput(1,3);   % Number of nodes per plane
    clear PlaneInput;
end

for nF = 9   %Number of phase
    Freq = 4;%Number of frequency
    tic
    display(['Phase ' num2str(nF)])
    Kmu=zeros([sum(NumP) 1]); % Allocate matrix with Mu and Kapa kernels
    Kka=zeros([sum(NumP) 1]);
    In=1;   %Starting index to sort multiple planes
    
    for ii=1:3:size(Planes,1)
        dKmu=zeros([1 1 NumP((ii+2)/3)]); % Allocate Mu and Kapa kernels
        dKka=zeros([1 1 NumP((ii+2)/3) 1]);
        
        for nPF=1:3
            PathF = [PathFor '/Source1/output/planes/planedisplacements.' ...
                num2str(ii+nPF-1)];
            PlaneF(nPF) = fopen(PathF);
            fread(PlaneF(nPF), 8,'double');  % Aux =
            fread(PlaneF(nPF), 1,'double');  % LengthY =
            fread(PlaneF(nPF), 1,'double');  % LengthX =
            fread(PlaneF(nPF), 1,'int');     % PointsDip =
            PointsStrike = fread(PlaneF(nPF), 1,'int');
            DeltaT = fread(PlaneF(nPF), 1,'double');
            Steps = fread(PlaneF(nPF), 1,'int');
            TimeSteps = round(t/(Steps*DeltaT))+1;
        end

        for nPA=1:3
            PathA = [PathAdj '/Phase' num2str(nF) '_Freq' num2str(Freq) ...
                '/output/planes/planedisplacements.' num2str(ii+nPA-1)];
            PlaneA(nPA) = fopen(PathA);
            fread(PlaneA(nPA), 8,'double');  % Aux =
            fread(PlaneA(nPA), 1,'double');  % LengthY =
            fread(PlaneA(nPA), 1,'double');  % LengthX =
            fread(PlaneA(nPA), 1,'int');     % PointsDip =
            fread(PlaneA(nPA), 1,'int');     % PointsStrike =
            fread(PlaneA(nPA), 1,'double');  % DeltaT =
            fread(PlaneA(nPA), 1,'int');     % Steps =
            fseek(PlaneA(nPA),24*TimeSteps*PointsStrike,0); % Go to the end of the file
        end
        
        % Allocate auxiliary deformation components
        U=zeros([3 3 PointsStrike]); V=U; 
        
        for nM=1:TimeSteps
            if mod(nM,20)==0
                disp(['Plane ' num2str(ii) '-' num2str(ii+2) ' of '   ...
                    num2str(size(Planes,1)) ', T = ' num2str(nM*Steps*DeltaT)])
            end
            
            for nPl=1:3
                fseek(PlaneA(nPl),-24*PointsStrike,0);
            end

            for nPl=1:3
                FP=fread(PlaneF(nPl), 3*PointsStrike,'double');
                AP=fread(PlaneA(nPl), 3*PointsStrike,'double');

                U(1,nPl,:) = FP(1:3:end);
                U(2,nPl,:) = FP(2:3:end);
                U(3,nPl,:) = FP(3:3:end);

                V(1,nPl,:) = AP(1:3:end);
                V(2,nPl,:) = AP(2:3:end);
                V(3,nPl,:) = AP(3:3:end);
            end
            
            % Compute Kmu and Kka at each time step
            dKka=dKka+(U(1,1,:)+U(2,2,:)+U(3,3,:)).*...
                      (V(1,1,:)+V(2,2,:)+V(3,3,:));
           
            dKmu=dKmu+0.5*((U(1,2,:)+U(2,1,:)).*(V(2,1,:)+V(1,2,:))+...
                           (U(1,3,:)+U(3,1,:)).*(V(3,1,:)+V(1,3,:))+...
                           (U(2,3,:)+U(3,2,:)).*(V(3,2,:)+V(2,3,:)));
            
            for nPl=1:3
                fseek(PlaneA(nPl),-24*PointsStrike,0);
            end
        end

        for nPl=1:3
            fclose(PlaneF(nPl));
            fclose(PlaneA(nPl));
        end
        
        % Sort Kmu and Kka en dK
        Kmu(In:In+NumP((ii+2)/3)-1)=squeeze(dKmu);
        Kka(In:In+NumP((ii+2)/3)-1)=squeeze(dKka);
        In=In+NumP((ii+2)/3);
        clear PlaneA PlaneF AP FP;
    end
      
    Kmu=(Steps*DeltaT)*Kmu*Factor/1e20; % Correct for forward simulation Factor, and integral and save
    Kka=(Steps*DeltaT)*Kka*Factor/1e20; % adjoint simulation 1e20 factor and integral DeltaT factor 
    save ([PathAdj 'KernelPhase' num2str(nF) 'Freq' num2str(Freq) '.mat'],'Kmu') ;
    save ([PathAdj 'KernelPhase' num2str(nF) 'Freq' num2str(Freq) '.mat'],'Kka','-append') ;
    
    clear Kmu Kka dKmu dKka U V;
       
end


