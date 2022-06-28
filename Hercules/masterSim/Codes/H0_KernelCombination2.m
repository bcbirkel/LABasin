% Kernel combination
% Written by Alan Juarez (USC) Feb, 2018.
clear; close all; clc;

% Model parameters
load basin.mat;

nP = 1:10;   % Phase number
nF = 1:3;   % Phase freq

% Load combination factors
load('InterFactFrec.mat');

% Load kernels
for i = 1:length(nP)
    for j=  1:length(nF)
        Phase(i,1).Freq(j,1).Kernel = load(['KernelP' num2str(nP(i)) 'F' num2str(nF(j)) '.mat']);
    end
end

% Define Kf the objective function
Kf = 0*vs(1:7.5e3/dxz+1,2+10e3/dxz:2+30e3/dxz);
Kf(abs(vs(1:7.5e3/dxz+1,2+10e3/dxz:2+30e3/dxz)-1420) < 100) = -1;
nz = size(Kf,1);
nx = size(Kf,2);
Kf = reshape(Kf,nx*nz,1);
Kf = [Kf; 0*Kf];

G = zeros([length(Kf) length(nP)*length(nF)]);

k = 1; 
for i = 1:length(nP)
    for j=  1:length(nF)
        [Dx,Dz] = gradient(Phase(i,1).Freq(j,1).Kernel.Kernel);
        G1 = (Dx+Dz)/dxz^2;
        G(:,k) = [reshape(C(i,i,j)*Phase(i,1).Freq(j,1).Kernel.Kernel,nx*nz,1) ...
            ; reshape(G1,nx*nz,1)];
        k = k+1;
    end
end

Alpha = (G'*G)\(G'*Kf);

Alpha = Alpha/max(abs(Alpha));


subplot(2,2,4)
for j=  1:length(nF)
    plot(Alpha(j:3:end),'-..','linewidth',1,'markersize',20), 
    hold on
end
axis 'square', xlim([0 length(nP)+1]), ylim([-1.1 1.1])
xlabel('Phase number'), ylabel('\alpha'), legend('2 Hz','3 Hz','4 Hz')
title 'Combination coefficient'


% Combine kernels
BASE_C = 0*BASE(i,:);
Ks = 0*Phase(1,1).Freq(1,1).Kernel.Kernel;

k = 1;
for i = 1:length(nP)
    for j=  1:length(nF)
            Ks = Ks+Alpha(k)*C(i,i,j)*Phase(i,1).Freq(j,1).Kernel.Kernel;
            BASE_C = BASE_C+Alpha(k)*BASE(i,:);
            k = k+1;
    end
end

%% Now you can start plots
figure;
imagesc([0:size(Ks,2)+1]*dxz/1e3-10,[0:size(Ks,1)+1]*dxz/1e3, Ks); hold on
plot(-6,2,'p','markerfacecolor',[1 0.84 0],...
    'markeredgecolor',[0 0 0],'MarkerSize',12,'color','k');
plot(6,0.1,'v','markerfacecolor',[0.47 0.67 0.19],...
    'markeredgecolor',[0 0 0],'MarkerSize',10,'color','k');
title(['Assembled  Kernel'])
ylabel('Depth (km)'); xlabel('Distance (km)'); axis equal tight;
CM = flip(loadcmap('BlueWhiteOrangeRed.c3g')); colormap(CM);
c1 = colorbar('Location','eastoutside'); title(c1,'K_v_s(s)');
Lim1 = 3*std(Ks(:));
caxis(Lim1*[-1 1]);

figure;
imagesc([0:size(Ks,2)+1]*dxz/1e3-10,[0:size(Ks,1)+1]*dxz/1e3, reshape(Kf(1:nz*nx),nz,nx)); hold on
title(['Target  Kernel'])
ylabel('Depth (km)'); xlabel('Distance (km)'); axis equal tight;
CM = flip(loadcmap('BlueWhiteOrangeRed.c3g')); colormap(CM);
c1 = colorbar('Location','eastoutside'); title(c1,'K_v_s(s)');
caxis([-1 1]);

%% Plot Seismogram Decomposition %%%
dt=5*0.004;     %Simulations dt
t=0:dt:12;  %Same as simulations

figure;
subplot(2,2,1)
plot(t,sum(BASE,1)/max(sum(BASE,1)),'k','linewidth',2), hold on

subplot(2,2,1)
plot(t,BASE_C/max(BASE_C),'r','linewidth',1.5),
legend('Data','Localized','location','southeast')
xlim([0 t(end)]), ylim(1.5*[-1 1]),
title('Localization')
xlabel('Time [s]') , ylabel('Displacement [m]')
