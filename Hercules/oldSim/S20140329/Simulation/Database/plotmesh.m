%
%   Script to plot the output mesh from Hercules:
%
%%

clear all; close all, clc
datum=1;

%%
FileID = fopen(['/home/dude/alanjuar/Documents/S20100707/Codes/AuxFiles/XYtoLonLat.txt']);

Zone = fscanf(FileID, '%f', 1); Az = fscanf(FileID, '%f', 1);
X0_lat = fscanf(FileID, '%f', 1); Y0_lon = fscanf(FileID, '%f', 1);
XLim(1) = fscanf(FileID, '%f', 1); XLim(2) = fscanf(FileID, '%f', 1);
YLim(1) = fscanf(FileID, '%f', 1); YLim(2) = fscanf(FileID, '%f', 1);
ZLim(1) = fscanf(FileID, '%f', 1); ZLim(2) = fscanf(FileID, '%f', 1);
fclose(FileID);

xi(1)=YLim(1); yi(1)=XLim(1);
xi(2)=YLim(1); yi(2)=XLim(2);
xi(3)=YLim(2); yi(3)=XLim(2);
xi(4)=YLim(2); yi(4)=XLim(1);
xi(5)=YLim(1); yi(5)=XLim(1);

figure;
plot3(xi*0.001,yi*0.001,yi*0+0,'k','LineWidth',1)
hold on


%%
fac=[ 1 2 3 4; 1 2 6 5; 2 3 7 6; 4 3 7 8; 5 6 7 8; 1 4 8 5]  ;


totElements=0;


for iFile=0:7
    iFile
    fp=fopen(['mesh_' num2str(iFile) '.mat'],'r');
    
    numElements=fread(fp,1,'int');
    totElements=totElements+numElements;
    aux=fread(fp,28*numElements,'double');
    aux=0.001*aux;
    for iE=1:2:numElements;
        
        i1=28*(iE-1);
        i2=28*(iE);
        elementS=[aux(i1+9:i1+16) aux(i1+1:i1+8) aux(i1+17:i1+24)];
        elementS(3,1:2)=elementS(8,1:2);
        elementS(4,1:2)=elementS(7,1:2);
        elementS(7:8,1:2)=elementS(3:4,1:2);
        
        element=elementS';
        vp=aux(i1+25);
        
        patch('faces',fac  ,'vertices',element','FaceColor',...
            [1 1 1],'EdgeColor',[.1 .1 .1 ],'linewidth',1.5)%,'FaceAlpha',.1)
        if(mod(iE,10000)==0)
            drawnow
            
            axis([0 200 0 200 -50 50])
        end
    end
    fclose(fp);
    drawnow
end
ylabel('Dist to North [Km]')
xlabel('Dist to East [Km]')
zlabel('Depth [Km]')
axis equal
set(gca,'zdir','reverse')
% axis equal
% view(3)

