% Process CVM426 and adapt it for Hercules
% Written by Alan Juarez (UNAM-USC) Dec, 2016.
clear; close all; clc;

ZMap=0;

FileID = fopen('cvms426_final_model.txt');
Header=fgetl(FileID);

N = 1;  Depth(N)=ZMap;
while ~feof(FileID) 
    X(N,1)=fscanf(FileID, '%f', 1);
    if ~isempty(X)
        Y(N,1)=fscanf(FileID, '%f', 1);
        Z(N,1)=fscanf(FileID, '%f', 1);
        fscanf(FileID, '%f', 4);
        Depth(N,1)=fscanf(FileID, '%f', 1);
        fscanf(FileID, '%f', 3);
        
        if Depth(N)-ZMap > 1
            break;
        end
        if abs(Depth(N)-ZMap)<1
            N = N + 1;
        end
        
    end
end
fclose(FileID);

x=unique(X);
y=unique(Y);

lx=length(x);
ly=length(y);
Z_Depth = 0:500:49500;
lz=length(Z_Depth);

FileID = fopen('cvms426_final_model.txt');
fgetl(FileID);

for iz=1:lz 
    for iy=1:ly
        for ix=1:lx
            fscanf(FileID, '%f', 3);
            fscanf(FileID, '%f', 4);
            
            Dh(ix,iy,iz)=fscanf(FileID, '%f', 1);
            Vp(ix,iy,iz)=fscanf(FileID, '%f', 1);
            Vs(ix,iy,iz)=fscanf(FileID, '%f', 1);
            Rho(ix,iy,iz)=fscanf(FileID, '%f', 1);
        end
    end
%     surf(Vs(:,:,10)), shading interp, view(2), axis equal;
end
fclose(FileID);

DFH=4; %Decimation factor in horizontal direction
DFV=2; %Decimation factor in vertical direction
File='vel_profiles_unidad1.fun';
FileW = fopen(File,'w');
fprintf(FileW,'%i\t',ly*lx/(DFH*DFH));

for iy=1:DFH:ly
    for ix=1:DFH:lx
        fprintf(FileW,'\n');
        fprintf(FileW,'%i\n',length(Dh(ix,iy,1:DFV:end)));
        fprintf(FileW,'%6.2f\t',Dh(ix,iy,1:DFV:end));
        fprintf(FileW,'\n');
        fprintf(FileW,'%6.2f\t',Vp(ix,iy,1:DFV:end));
        fprintf(FileW,'\n');
        fprintf(FileW,'%6.2f\t',Vs(ix,iy,1:DFV:end));
        fprintf(FileW,'\n');
        fprintf(FileW,'%6.2f\t',Rho(ix,iy,1:DFV:end));
    end
end
fclose(FileID);

P=0:lx*ly/(DFH*DFH)-1;
File='vel_profiles_unidad1.surf';
fileID = fopen(File,'w');
fprintf(fileID,'%i\t',[ly/DFH lx/DFH]);
fprintf(fileID,'\n');
fprintf(fileID,'%6.2f\t',500*(y(1:DFH:end)-1));
fprintf(fileID,'\n');
fprintf(fileID,'%6.2f\t',500*(x(1:DFH:end)-1));
fprintf(fileID,'\n');
fprintf(fileID,'%i\t',P);
fclose(fileID); 