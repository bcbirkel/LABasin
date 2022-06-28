function dat=display_globec_grid;

linear_scale=0;  % 0 = log scale, 1 = linear scale

[filename, filepath]=uigetfile('*.mat','Select An Input Data File');

if isempty(filename)
   return
end

fname=[filepath filename];

load(fname)

clr_max=length(colormap);

X=data.out.krig.Xg;
Y=data.out.krig.Yg;
x=data.out.krig.gx;         % GLOBEC grids in x direction
y=data.out.krig.gy;         % GLOBEC grids in y direction
if linear_scale == 1
   V=exp(data.out.krig.Vg);   % convert log scale to linear scale
   E=exp(data.out.krig.Eg);   % 
   var=exp(data.out.krig.gv)-1;
else
   V=data.out.krig.Vg;      % krig outputs are ln(z+1)
   E=data.out.krig.Eg;      % krig outputs are ln(z+1)
   var=data.out.krig.gv;   
end
n=length(x);

%% Plot kriging results on regular grids
figure(1)
pcolor(X,Y,V);
shading interp;colorbar
hold on
[c,H]=contour(X,Y,V,'k');
hh=clabel(c,H);
hold off
xlabel('LONGITUDE','fontsize',16,'fontweight','bold')  
ylabel('LATITUDE','fontsize',16,'fontweight','bold')
title('KRIGING MAP ON REGULAR GRIDS','fontsize',16,'fontweight','bold');

figure(2)
pcolor(X,Y,E);
hold on
shading interp;colorbar
[c,H]=contour(X,Y,E,'k');xlabel('LONGITUDE','fontsize',16,'fontweight','bold')  
hh=clabel(c,H);
hold off
ylabel('LATITUDE','fontsize',16,'fontweight','bold')
title('KRIGING VARIANCE ON REGULAR GRIDS','fontsize',16,'fontweight','bold');
colorbar

%% Plot kriging resluts on the GLOBEC grids
clr=max(0,floor(clr_max*(var-min(var))/(max(var)-min(var))))+1;
indx_nan=find(isnan(var) == 1);
clr=min(clr,clr_max);
cmap=colormap;
figure(3)
i0=1;
for i=1:n
  if isempty(indx_nan) | min(abs(indx_nan-i)) ~= 0
     plot(x(i),y(i),'.','color',cmap(clr(i),:),'markersize',14);
     if i == i0;hold on,end
  else
     i0=i0+1;
  end
end
xlabel('LONGITUDE','fontsize',16,'fontweight','bold')  
ylabel('LATITUDE','fontsize',16,'fontweight','bold')
title('KRIGING RESULTS ON THE GLOBEC GRIDS','fontsize',16,'fontweight','bold');
pos = get(gca,'Position'); 
stripe = 0.075; edge = 0.02; 
[az,el] = view;
if all([az,el]==[0 90]), space = 0.05; else space = .1; end
set(gca,'Position',[pos(1) pos(2) pos(3)*(1-stripe-edge-space) pos(4)])
rect = [pos(1)+(1-stripe-edge)*pos(3) pos(2) stripe*pos(3) pos(4)];
ax= axes('Position', rect);
image([0 1],[min(var) max(var)],(1:clr_max)','Tag','TMW_COLORBAR','deletefcn','colorbar(''delete'')'); 
set(ax,'Ydir','normal')
set(ax,'YAxisLocation','right')
set(ax,'xtick',[])
dat=data;
  