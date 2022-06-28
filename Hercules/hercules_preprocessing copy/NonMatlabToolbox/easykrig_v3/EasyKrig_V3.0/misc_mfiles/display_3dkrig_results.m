% display 3D kriging results

addpath C:\projects\EasyKrig\EasyKrig_V3.0_Beta\bin C:\projects\EasyKrig\EasyKrig_V3.0_Beta\general

clear

clf
[filename, filepath]=uigetfile('*.mat','Select An Input Data File');

if isempty(filename)
   return
end

fname=[filepath filename];

load(fname)

X=data.out.krig.Xg;
Y=data.out.krig.Yg;
Z=data.out.krig.Zg;
V=data.out.krig.Vg;
E=data.out.krig.Eg;

xm=squeeze(mean(X(:,:,1)))';
ym=squeeze(Y(:,1,1));
zm=squeeze(Z(1,1,:));

slice_index_x=9;
slice_index_y=1;
slice_index_z=19;

ntick=4;

figure(1)
slice(X,Y,Z,V,[xm(slice_index_x)],[ym(slice_index_y)],[zm(slice_index_z)])
xlabel('LONGITUDE (deg)','fontsize',16,'fontweight','bold')
ylabel('LATITUDE (deg)','fontsize',16,'fontweight','bold')
zlabel('DEPTH (m)','fontsize',16,'fontweight','bold')
title('TEMPERATURE (^oC)','fontsize',16,'fontweight','bold')
[xinc,xdits,yinc,ydits]=get_ninc(gca,ntick);
mapax(xinc,xdits,yinc,ydits,gca,1);		% x-axis label
mapax(xinc,xdits,yinc,ydits,gca,2);		% y-axis label
colorbar
set(gca,'zdir','reverse')
shading flat
axis([-22 -12.5 -25 -20 3500 6000])

figure(2)
slice(X,Y,Z,E,[xm(slice_index_x)],[ym(slice_index_y)],[zm(slice_index_z)])
xlabel('LONGITUDE (deg)','fontsize',16,'fontweight','bold')
ylabel('LATITUDE (deg)','fontsize',16,'fontweight','bold')
zlabel('DEPTH (m)','fontsize',16,'fontweight','bold')
title('KRIGING VARIANCE','fontsize',16,'fontweight','bold')
[xinc,xdits,yinc,ydits]=get_ninc(gca,ntick);
mapax(xinc,xdits,yinc,ydits,gca,1);		% x-axis label
mapax(xinc,xdits,yinc,ydits,gca,2);		% y-axis label
colorbar
set(gca,'zdir','reverse')
shading flat
axis([-22 -12.5 -25 -20 3500 6000])
