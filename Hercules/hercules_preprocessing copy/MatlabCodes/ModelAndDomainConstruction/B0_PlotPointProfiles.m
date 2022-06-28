% Process and plot CVM426
% Written by Alan Juarez (UNAM-USC) Dec, 2016.
clear; close all; clc;

% Load XYZ mesh data
Mesh = load(['/home/scec-00/alanjuar/EarthquakePhysics/Codes/AuxFiles/VelModelCoords.xyz']);
X = Mesh(:,1);
Y = Mesh(:,2);
clear Mesh;

% Load UCVM model
fileID = fopen('/home/scec-00/alanjuar/EarthquakePhysics/Codes/AuxFiles/ModelOut.out','r');
formatSpec = ['%f %f %f %f %f %s %f %f %f %s %f %f %f %s %f %f %f'];
A = textscan(fileID,formatSpec);
fclose(fileID);

% Assign properties to arrays
% Lon = A{1,1};
% Lat = A{1,2};
Z = A{1,3};
% Surf = A{1,4};
% vs30 = A{1,5};
% crustal = A{1,6};
% cr_vp = A{1,7};
% cr_vs = A{1,8};
% cr_rho = A{1,9};
% gtl = A{1,10};
% gtl_vp = A{1,11};
% gtl_vs = A{1,12};
% gtl_rho = A{1,13};
% cmb_algo = A{1,14};
Vp = A{1,15};
Vs = A{1,16};
Rho = A{1,17};
clear A;

% % % Estimate 30% velocity drop
% In = find(Vs <= 1200 & Z <= 1500);
% Vs(In) = 0.7*Vs(In);


% unique arrays
x = unique(X);
y = unique(Y);
z = unique(Z);
ly = length(y);
lx = length(x);
lz = length(z);

% Sort Data into 3D arrays
nP=1;
for ix=1:lx
    for iy=1:ly
        for iz = 1:lz
            Depth(ix,iy,iz) = Z(nP);
            Alpha(ix,iy,iz) = Vp(nP);
            Beta(ix,iy,iz) = Vs(nP);
            Density(ix,iy,iz) = Rho(nP);
            nP = nP+1;
        end
    end
end

% Open file for velprofiles.fun which contains vertical profiles with
% depths
ix = 101;
iy = 15;

Z1 = squeeze(Depth(ix,iy,:));
Vp1 = squeeze(Alpha(ix,iy,:));
Vs1 = squeeze(Beta(ix,iy,:));
Rho1 = squeeze(Density(ix,iy,:));


% Estimate velocity change
Vs2 = Vs1;
for n = 1:length(Z1)
    if Z1(n) < 600 && Vs1(n) < 1005
        Vs2(n) = (0.3*Z1(n)/500+0.7)*Vs1(n);
    end
end

figure;
plot(0.001*Vs1,0.001*Z1,0.001*Vp1,0.001*Z1,0.001*Rho1,0.001*Z1,'linewidth',2), hold on
plot(0.001*Vs2,0.001*Z1,'.','linewidth',2.5), axis 'ij'
ylim([0 5]), ylabel 'Depth(km)',xlabel 'Vp,Vs(km/s), \rho.(g/cm^{3})'
legend('Vs','Vp','\rho','0.7Vs'), grid minor

print('ProfileSurf','-dpng','-r800','-painters')

figure;
hold on
for ix=1:lx
    for iy=1:ly
        plot(0.001*squeeze(Beta(ix,iy,:)),0.001*Z1,'color',0.8*[1 1 1]),
        plot(0.001*squeeze(Alpha(ix,iy,:)),0.001*Z1,'color',0.6*[1 1 1]),
    end
end

m_Vs = 0.001*squeeze(mean(mean(Beta)));
plot(m_Vs,0.001*Z1,'linewidth',2.5),

m_Vp = 0.001*squeeze(mean(mean(Beta)));
plot(m_Vp,0.001*Z1,'linewidth',2.5), axis 'ij'
ylim([0 5]), ylabel 'Depth(km)',xlabel 'Vp,Vs(km/s)'
legend('Vs','Vp','\rho','0.7Vs'), grid minor
