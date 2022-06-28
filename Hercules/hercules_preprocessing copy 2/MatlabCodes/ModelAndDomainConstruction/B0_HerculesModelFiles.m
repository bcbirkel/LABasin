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
FileW = fopen('vel_profiles_unidad1.fun','w');
fprintf(FileW,'%i\t',ly*lx);

for ix=1:lx
    for iy=1:ly
        fprintf(FileW,'\n');
        fprintf(FileW,'%i\n',lz);
        fprintf(FileW,'%6.2f\t',Depth(ix,iy,:));
        fprintf(FileW,'\n');
        fprintf(FileW,'%6.2f\t',Alpha(ix,iy,:));
        fprintf(FileW,'\n');
        fprintf(FileW,'%6.2f\t',Beta(ix,iy,:));
        fprintf(FileW,'\n');
        fprintf(FileW,'%6.2f\t',Density(ix,iy,:));
    end
end
fclose(FileW);

% Open file for velprofiles.surf wich contains the locations of the
% profiles and the umber of profiles
P = 0:lx*ly-1;
fileX = fopen('vel_profiles_unidad1.surf','w');
fprintf(fileX,'%i\t',[lx ly]);
fprintf(fileX,'\n');
fprintf(fileX,'%6.2f\t',x);
fprintf(fileX,'\n');
fprintf(fileX,'%6.2f\t',y);
fprintf(fileX,'\n');
fprintf(fileX,'%i\t',P);
fclose(fileX); 
