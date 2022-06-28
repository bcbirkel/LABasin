% GRT_calculator.m
% code to compute generalized ray theory synthetics for a given source
% model and a given earth structure model
% usage:
% synth = GRT_calculator(model_params, source_params)
% code is far from robust, but will compute a series of rays for region
% scale model
%
% model_params contains the layer interface depths in the first row, layer
% shear velocities in the second row, and layer densities in the third row.
% example:
% model_params = [0,4,17,25;2.4,3.5,3.9,4.5;2.0,2.6,2.8,3.3]
%
% source_params contains various parameters relating to the earthquake
% source
% depth (km), distance (km), M0 (erg), rise_time (seconds), back_azimuth,
% strike, rake, dip
% example:
% source_params = [11, 63.2, 1.85*10^24, .25, 320.6, 323, 180, 90]
%
%Defining the variables
% z=[0,4,17,25]; %layer interface depths (km)
% h=11;           %Source depth (km)
% beta=[2.4,3.5,3.9,4.5];   %velocity of layer (km/s)
% r=63.8;    %lateral distance from source to receiver (km)
% rho=[2.0,2.6,2.8,3.3];%Density of layers (g/cm^3)
% %shear moduli of layers
% mu=zeros(length(beta));
% for i=1:length(beta)
%     mu(i)=rho(i)*beta(i)^2;
% end
% M0=1.85*10^24;   %Scalar seismic moment (erg)
% tau=.25;
% azi=320.6;
% strike=323;
% theta=(azi-strike)*pi/180;
% lambda=180*pi/180;
% delta=90*pi/180;
% dt=.01;

function v = GRT_calculator(model,source)
z=model(1,:);
beta=model(2,:);
rho=model(3,:);

h=source(1);
r=source(2);
M0=source(3);
tau=source(4);
azi=source(5);
strike=source(6);
lambda = source(7) * pi/180;
delta = source(8) * pi/180;
dt = 0.01;
mu = zeros(length(beta));
for i=1:length(beta);
    mu(i) = rho(i) * beta(i) * beta(i);
end
theta = (azi-strike) * pi/180;

%A parameters
A4=cos(2*theta)*cos(lambda)*sin(delta)-.5*sin(2*theta)*sin(lambda)*sin(2*delta);
A5=-sin(theta)*cos(lambda)*cos(delta)-cos(theta)*sin(lambda)*cos(2*delta);

%Scaling coefficient
scale=(2*M0*sqrt(2))/(4*pi^2*rho(2)*sqrt(r)*beta(2)^2)*10^-20;

%Source time history vector
tt=0:.01:2;
sth=tt.*exp(-tt./tau);
integsth=sum(sth.*.01);
sthnorm=sth./integsth;

%Defining search grid for complex t
%Direct Ray
prealmax=.6;
timedirect(1)=0;
nnpts=200;
dp=prealmax/nnpts;
preal(1)=0.05;
pimag(1)=0;
theta=0:pi/180:pi/2;
for i=2:nnpts
    for j=1:length(theta)
        prealdummy(j)=preal(i-1)+dp*cos(theta(j));
        pimagdummy(j)=pimag(i-1)+dp*sin(theta(j));
        pdummy(j)=prealdummy(j)+sqrt(-1)*pimagdummy(j);
        timedummy(j)=pdummy(j)*r+(h-z(2))*sqrt(1/beta(2)^2-pdummy(j)^2)+z(2)*sqrt(1/beta(1)^2-pdummy(j)^2);
        if j==1
            goldenapple=abs(imag(timedummy(j))); %goldenapple = 'golden apple' aka the value we want. Search to minimize goldenapple(j)
            l=j;
        else
            if abs(imag(timedummy(j)))<goldenapple
                goldenapple=abs(imag(timedummy(j)));
                l=j;
            end
        end
    end
    preal(i)=prealdummy(l);
    pimag(i)=pimagdummy(l);
    timedirect(i)=real(timedummy(l));
end
xx=0:dt:50;
hmm=sortrows([timedirect;preal;pimag]');
%Remove duplicates and backtracks
k=1;
for i=2:length(hmm(:,2))
    if hmm(i,2)>=hmm(i-1,2) && hmm(i,3)>=hmm(i-1,3)
        hmm(k,1)=hmm(i,1);
        hmm(k,2)=hmm(i,2);
        hmm(k,3)=hmm(i,3);
        k=k+1;
    end
end
i=2;
while i<=length(hmm(:,1))
    if hmm(i,1)<=hmm(i-1,1)
        hmm(i:length(hmm(:,1)),:)=[];
    else
        i=i+1;
    end    
end
pcomplex=hmm(:,2)+hmm(:,3)*sqrt(-1);
timedirect=hmm(:,1)';
pdirectcontour=interp1(timedirect,pcomplex,xx,'linear','extrap');
% figure(1)
% subplot(2,5,1), scatter(preal,pimag)
% hold on
% subplot(2,5,1), plot(real(pdirectcontour),imag(pdirectcontour))
% axis([0 1 0 .1])
clear hmm
%csvwrite('directrob.csv',pdirectcontour)

%Defining search grid for complex t
%Upmult1
timeupmult1(1)=0;
for i=2:nnpts
    for j=1:length(theta)
        prealdummy(j)=preal(i-1)+dp*cos(theta(j));
        pimagdummy(j)=pimag(i-1)+dp*sin(theta(j));
        pdummy(j)=prealdummy(j)+sqrt(-1)*pimagdummy(j);
        timeupmult1dummy(j)=pdummy(j)*r+(h-z(2))*sqrt(1/beta(2)^2-pdummy(j)^2)+3*z(2)*sqrt(1/beta(1)^2-pdummy(j)^2);
        if j==1
            goldenapple=abs(imag(timeupmult1dummy(j))); %goldenapple = 'golden apple' aka the value we want. Search to minimize GA(j)
            l=j;
        else
            if abs(imag(timeupmult1dummy(j)))<goldenapple
                goldenapple=abs(imag(timeupmult1dummy(j)));
                l=j;
            end
        end
    end
        preal(i)=prealdummy(l);
        pimag(i)=pimagdummy(l);
        timeupmult1(i)=real(timeupmult1dummy(l));
end
hmm=sortrows([timeupmult1;preal;pimag]');
k=1;
for i=2:length(hmm(:,2))
    if hmm(i,2)>=hmm(i-1,2) && hmm(i,3)>=hmm(i-1,3)
        hmm(k,1)=hmm(i,1);
        hmm(k,2)=hmm(i,2);
        hmm(k,3)=hmm(i,3);
        k=k+1;
    end
end
i=2;
while i<=length(hmm(:,1))
    if hmm(i,1)<=hmm(i-1,1)
        hmm(i:length(hmm(:,1)),:)=[];
    else
        i=i+1;
    end    
end
pcomplex=hmm(:,2)+hmm(:,3)*sqrt(-1);
timeupmult1=hmm(:,1)';
pupmult1contour=interp1(timeupmult1,pcomplex,xx,'linear','extrap');
% subplot(2,5,2), scatter(preal,pimag)
% hold on
% subplot(2,5,2), plot(real(pupmult1contour),imag(pupmult1contour))
% axis([0 1 0 .2])
clear hmm

%Defining search grid for complex t
%upmult 2
timeupmult2(1)=0;
for i=2:nnpts
    for j=1:length(theta)
        prealdummy(j)=preal(i-1)+dp*cos(theta(j));
        pimagdummy(j)=pimag(i-1)+dp*sin(theta(j));
        pdummy(j)=prealdummy(j)+sqrt(-1)*pimagdummy(j);
        timeupmult2dummy(j)=pdummy(j)*r+(h-z(2))*sqrt(1/beta(2)^2-pdummy(j)^2)+5*z(2)*sqrt(1/beta(1)^2-pdummy(j)^2);
        if j==1
            goldenapple=abs(imag(timeupmult2dummy(j))); %goldenapple = 'golden apple' aka the value we want. Search to minimize GA(j)
            l=j;
        else
            if abs(imag(timeupmult2dummy(j)))<goldenapple
                goldenapple=abs(imag(timeupmult2dummy(j)));
                l=j;
            end
        end
    end
        preal(i)=prealdummy(l);
        pimag(i)=pimagdummy(l);
        timeupmult2(i)=real(timeupmult2dummy(l));
end
hmm=sortrows([timeupmult2;preal;pimag]');
k=1;
for i=2:length(hmm(:,2))
    if hmm(i,2)>=hmm(i-1,2) && hmm(i,3)>=hmm(i-1,3) && hmm(i,1)>=hmm(i-1,3)
        hmm(k,1)=hmm(i,1);
        hmm(k,2)=hmm(i,2);
        hmm(k,3)=hmm(i,3);
        k=k+1;
    end
end
i=2;
while i<=length(hmm(:,1))
    if hmm(i,1)<=hmm(i-1,1)
        hmm(i:length(hmm(:,1)),:)=[];
    else
        i=i+1;
    end    
end
pcomplex=hmm(:,2)+hmm(:,3)*sqrt(-1);
timeupmult2=hmm(:,1)';
pupmult2contour=interp1(timeupmult2,pcomplex,xx,'linear','extrap');
% subplot(2,5,3), scatter(preal,pimag)
% hold on
% subplot(2,5,3), plot(real(pupmult2contour),imag(pupmult2contour))
% axis([0 1 0 .2])
clear hmm
%Defining search grid for complex t
%up mult 3
timeupmult3(1)=0;
for i=2:nnpts
    for j=1:length(theta)
        prealdummy(j)=preal(i-1)+dp*cos(theta(j));
        pimagdummy(j)=pimag(i-1)+dp*sin(theta(j));
        pdummy(j)=prealdummy(j)+sqrt(-1)*pimagdummy(j);
        timeupmult3dummy(j)=pdummy(j)*r+(h-z(2))*sqrt(1/beta(2)^2-pdummy(j)^2)+7*z(2)*sqrt(1/beta(1)^2-pdummy(j)^2);
        if j==1
            goldenapple=abs(imag(timeupmult3dummy(j))); %goldenapple = 'golden apple' aka the value we want. Search to minimize GA(j)
            l=j;
        else
            if abs(imag(timeupmult3dummy(j)))<goldenapple
                goldenapple=abs(imag(timeupmult3dummy(j)));
                l=j;
            end
        end
    end
        preal(i)=prealdummy(l);
        pimag(i)=pimagdummy(l);
        timeupmult3(i)=real(timeupmult3dummy(l));
end
hmm=sortrows([timeupmult3;preal;pimag]');
k=1;
for i=2:length(hmm(:,2))
    if hmm(i,2)>=hmm(i-1,2) && hmm(i,3)>=hmm(i-1,3)
        hmm(k,1)=hmm(i,1);
        hmm(k,2)=hmm(i,2);
        hmm(k,3)=hmm(i,3);
        k=k+1;
    end
end
i=2;
while i<=length(hmm(:,1))
    if hmm(i,1)<=hmm(i-1,1)
        hmm(i:length(hmm(:,1)),:)=[];
    else
        i=i+1;
    end    
end
pcomplex=hmm(:,2)+hmm(:,3)*sqrt(-1);
timeupmult3=hmm(:,1)';
pupmult3contour=interp1(timeupmult3,pcomplex,xx,'linear','extrap');
% subplot(2,5,4), scatter(preal,pimag)
% hold on
% subplot(2,5,4), plot(real(pupmult3contour),imag(pupmult3contour))
% axis([0 1 0 .25])

%Defining search grid for complex t
%downref1
timedownref1(1)=0;
for i=2:nnpts
    for j=1:length(theta)
        prealdummy(j)=preal(i-1)+dp*cos(theta(j));
        pimagdummy(j)=pimag(i-1)+dp*sin(theta(j));
        pdummy(j)=prealdummy(j)+sqrt(-1)*pimagdummy(j);
        timedownref1dummy(j)=pdummy(j)*r+(z(3)-h+z(3)-z(2))*sqrt(1/beta(2)^2-pdummy(j)^2)+z(2)*sqrt(1/beta(1)^2-pdummy(j)^2);
        if j==1
            goldenapple=abs(imag(timedownref1dummy(j))); %goldenapple = 'golden apple' aka the value we want. Search to minimize GA(j)
            l=j;
        else
            if abs(imag(timedownref1dummy(j)))<goldenapple
                goldenapple=abs(imag(timedownref1dummy(j)));
                l=j;
            end
        end
    end
        preal(i)=prealdummy(l);
        pimag(i)=pimagdummy(l);
        timedownref1(i)=real(timedownref1dummy(l));
end
hmm=sortrows([timedownref1;preal;pimag]');
k=1;
for i=2:length(hmm(:,2))
    if hmm(i,2)>=max(hmm(1:i-1,2)) && hmm(i,3)>=max(hmm(1:i-1,3)) && hmm(i,1)>=max(hmm(1:i-1,1))
        hmm(k,1)=hmm(i,1);
        hmm(k,2)=hmm(i,2);
        hmm(k,3)=hmm(i,3);
        k=k+1;
    end
end
i=2;
while i<=length(hmm(:,1))
    if hmm(i,1)<=hmm(i-1,1)
        hmm(i:length(hmm(:,1)),:)=[];
    else
        i=i+1;
    end    
end
pcomplex=hmm(:,2)+hmm(:,3)*sqrt(-1);
timedownref1=hmm(:,1)';
pdownref1contour=interp1(timedownref1,pcomplex,xx,'linear','extrap');
% subplot(2,5,5), scatter(preal,pimag)
% hold on
% subplot(2,5,5), plot(real(pdownref1contour),imag(pdownref1contour))
% axis([0 1 0 .2])
%csvwrite('downref1contourrob.csv',pdownref1contour);
%csvwrite('downref1rawrob.csv',pcomplex);

%Downref 2
timedownref2(1)=0;
for i=2:nnpts
    for j=1:length(theta)
        prealdummy(j)=preal(i-1)+dp*cos(theta(j));
        pimagdummy(j)=pimag(i-1)+dp*sin(theta(j));
        pdummy(j)=prealdummy(j)+sqrt(-1)*pimagdummy(j);
        timedownref2dummy(j)=pdummy(j)*r+(z(3)-h+(z(3)-z(2)))*sqrt(1/beta(2)^2-pdummy(j)^2)+2*(z(4)-z(3))*sqrt(1/beta(3)^2-pdummy(j)^2)+z(2)*sqrt(1/beta(1)^2-pdummy(j)^2);
        if j==1
            goldenapple=abs(imag(timedownref2dummy(j))); %goldenapple = 'golden apple' aka the value we want. Search to minimize GA(j)
            l=j;
        else
            if abs(imag(timedownref2dummy(j)))<goldenapple
                goldenapple=abs(imag(timedownref2dummy(j)));
                l=j;
            end
        end
    end
        preal(i)=prealdummy(l);
        pimag(i)=pimagdummy(l);
        timedownref2(i)=real(timedownref2dummy(l));
end
%xx=min(timedownref2):dt:max(timedownref2);
hmm=sortrows([timedownref2;preal;pimag]');
k=1;
for i=2:length(hmm(:,2))
    if hmm(i,2)>=max(hmm(1:i-1,2)) && hmm(i,3)>=max(hmm(1:i-1,3)) && hmm(i,1)>=max(hmm(1:i-1,1))
        hmm(k,1)=hmm(i,1);
        hmm(k,2)=hmm(i,2);
        hmm(k,3)=hmm(i,3);
        k=k+1;
    end
end
i=2;
while i<=length(hmm(:,1))
    if hmm(i,1)<=hmm(i-1,1)
        hmm(i:length(hmm(:,1)),:)=[];
    else
        i=i+1;
    end    
end
pcomplex=hmm(:,2)+hmm(:,3)*sqrt(-1);
timedownref2=hmm(:,1)';
pdownref2contour=interp1(timedownref2,pcomplex,xx,'linear','extrap');
% subplot(2,5,6), scatter(preal,pimag)
% hold on
% subplot(2,5,6), plot(real(pdownref2contour),imag(pdownref2contour))
% %axis([0 1 0 .2])

%up to midcrust reflection
timeupdown1(1)=0;
for i=2:nnpts
    for j=1:length(theta)
        prealdummy(j)=preal(i-1)+dp*cos(theta(j));
        pimagdummy(j)=pimag(i-1)+dp*sin(theta(j));
        pdummy(j)=prealdummy(j)+sqrt(-1)*pimagdummy(j);
        timeupdown1dummy(j)=pdummy(j)*r+(h-z(2)+2*(z(3)-z(2)))*sqrt(1/beta(2)^2-pdummy(j)^2)+z(2)*sqrt(1/beta(1)^2-pdummy(j)^2);
        if j==1
            goldenapple=abs(imag(timeupdown1dummy(j))); %goldenapple = 'golden apple' aka the value we want. Search to minimize GA(j)
            l=j;
        else
            if abs(imag(timeupdown1dummy(j)))<goldenapple
                goldenapple=abs(imag(timeupdown1dummy(j)));
                l=j;
            end
        end
    end
        preal(i)=prealdummy(l);
        pimag(i)=pimagdummy(l);
        timeupdown1(i)=real(timeupdown1dummy(l));
end
hmm=sortrows([timeupdown1;preal;pimag]');
k=1;
for i=2:length(hmm(:,2))
    if hmm(i,2)>=max(hmm(1:i-1,2)) && hmm(i,3)>=max(hmm(1:i-1,3)) && hmm(i,1)>=max(hmm(1:i-1,1))
        hmm(k,1)=hmm(i,1);
        hmm(k,2)=hmm(i,2);
        hmm(k,3)=hmm(i,3);
        k=k+1;
    end
end
i=2;
while i<=length(hmm(:,1))
    if hmm(i,1)<=hmm(i-1,1)
        hmm(i:length(hmm(:,1)),:)=[];
    else
        i=i+1;
    end    
end
pcomplex=hmm(:,2)+hmm(:,3)*sqrt(-1);
timeupdown1=hmm(:,1)';
pupdown1contour=interp1(timeupdown1,pcomplex,xx,'linear','extrap');
% subplot(2,5,7), scatter(preal,pimag)
% hold on
% subplot(2,5,7), plot(real(pupdown1contour),imag(pupdown1contour))
%axis([0 1 0 .2])

%up to moho reflection
timeupdown2(1)=0;
for i=2:nnpts
    for j=1:length(theta)
        prealdummy(j)=preal(i-1)+dp*cos(theta(j));
        pimagdummy(j)=pimag(i-1)+dp*sin(theta(j));
        pdummy(j)=prealdummy(j)+sqrt(-1)*pimagdummy(j);
        timeupdown2dummy(j)=pdummy(j)*r+(h-z(2)+(z(3)-z(2))*2)*sqrt(1/beta(2)^2-pdummy(j)^2)+(2*(z(4)-z(3)))*sqrt(1/beta(3)^2-pdummy(j)^2)+z(2)*sqrt(1/beta(1)^2-pdummy(j)^2);
        if j==1
            goldenapple=abs(imag(timeupdown2dummy(j))); %goldenapple = 'golden apple' aka the value we want. Search to minimize GA(j)
            l=j;
        else
            if abs(imag(timeupdown2dummy(j)))<goldenapple
                goldenapple=abs(imag(timeupdown2dummy(j)));
                l=j;
            end
        end
    end
        preal(i)=prealdummy(l);
        pimag(i)=pimagdummy(l);
        timeupdown2(i)=real(timeupdown2dummy(l));
end
hmm=sortrows([timeupdown2;preal;pimag]');
k=1;
for i=2:length(hmm(:,2))
    if hmm(i,2)>=max(hmm(1:i-1,2)) && hmm(i,3)>=max(hmm(1:i-1,3)) && hmm(i,1)>=max(hmm(1:i-1,1))
        hmm(k,1)=hmm(i,1);
        hmm(k,2)=hmm(i,2);
        hmm(k,3)=hmm(i,3);
        k=k+1;
    end
end
i=2;
while i<=length(hmm(:,1))
    if hmm(i,1)<=hmm(i-1,1)
        hmm(i:length(hmm(:,1)),:)=[];
    else
        i=i+1;
    end    
end
pcomplex=hmm(:,2)+hmm(:,3)*sqrt(-1);
timeupdown2=hmm(:,1)';
pupdown2contour=interp1(timeupdown2,pcomplex,xx,'linear','extrap');
% subplot(2,5,8), scatter(preal,pimag)
% hold on
% subplot(2,5,8), plot(real(pupdown2contour),imag(pupdown2contour))
% %axis([0 1 0 .2])

%computing rt
A1=pi*pi/16;
sdt=2.0/sqrt(dt);
a=0;

for I=1:length(xx)
   b=sqrt((I-1) + A1);
   rt(I)=sdt * (b - a);
   a=b;
end

%preallocating eta functions
eta1=zeros(length(pdirectcontour));
eta2=zeros(length(pdirectcontour));
eta3=zeros(length(pdirectcontour));
eta4=zeros(length(pdirectcontour));
% 

dp(1)=0;
Jdirect1(1)=0;
Jdirect2(1)=0;
for i=2:length(pdirectcontour)
    %Direct
    eta1(i)=sqrt(1/beta(1)^2-pdirectcontour(i)^2);
    eta2(i)=sqrt(1/beta(2)^2-pdirectcontour(i)^2);
    eta3(i)=sqrt(1/beta(3)^2-pdirectcontour(i)^2);
    eta4(i)=sqrt(1/beta(4)^2-pdirectcontour(i)^2);
    transmission(i)=(2*(mu(2)*eta2(i))/(mu(2)*eta2(i)+mu(1)*eta1(i)));
    PI1(i)=transmission(i);
    dp(i)=(pdirectcontour(i)-pdirectcontour(i-1))/dt;
    Jdirect1(i)=imag(PI1(i)*pdirectcontour(i)^1.5*dp(i)/eta2(i));
    Jdirect2(i)=imag(PI1(i)*pdirectcontour(i)^0.5*dp(i));
end

for i=2:length(pupmult1contour)
    %up mult 1
     eta1(i)=sqrt(1/beta(1)^2-pupmult1contour(i)^2);
     eta2(i)=sqrt(1/beta(2)^2-pupmult1contour(i)^2);
     eta3(i)=sqrt(1/beta(3)^2-pupmult1contour(i)^2);
     eta4(i)=sqrt(1/beta(4)^2-pupmult1contour(i)^2);
     transmission(i)=(2*(mu(2)*eta2(i))/(mu(2)*eta2(i)+mu(1)*eta1(i)));
     reflection(i)=((mu(1)*eta1(i)-mu(2)*eta2(i))/(mu(2)*eta2(i)+mu(1)*eta1(i)))^1;
     PI2(i)=reflection(i)*transmission(i);
     dp(i)=(pupmult1contour(i)-pupmult1contour(i-1))/dt;
     Jupmult11(i)=imag(PI2(i)*pupmult1contour(i)^1.5*dp(i)/eta2(i));
     Jupmult12(i)=imag(PI2(i)*pupmult1contour(i)^0.5*dp(i));
end %for

for i=2:length(pupmult2contour)
    %up mult 2
     eta1(i)=sqrt(1/beta(1)^2-pupmult2contour(i)^2);
     eta2(i)=sqrt(1/beta(2)^2-pupmult2contour(i)^2);
     eta3(i)=sqrt(1/beta(3)^2-pupmult2contour(i)^2);
     eta4(i)=sqrt(1/beta(4)^2-pupmult2contour(i)^2);
     transmission(i)=(2*(mu(2)*eta2(i))/(mu(2)*eta2(i)+mu(1)*eta1(i)));
     reflection(i)=((mu(1)*eta1(i)-mu(2)*eta2(i))/(mu(2)*eta2(i)+mu(1)*eta1(i)))^2;
     PI3(i)=reflection(i)*transmission(i);
     dp(i)=(pupmult2contour(i)-pupmult2contour(i-1))/dt;
     Jupmult21(i)=imag(PI3(i)*pupmult2contour(i)^1.5*dp(i)/eta2(i));
     Jupmult22(i)=imag(PI3(i)*pupmult2contour(i)^0.5*dp(i));
end %for

for i=2:length(pupmult3contour)
    %up mult 3
     eta1(i)=sqrt(1/beta(1)^2-pupmult3contour(i)^2);
     eta2(i)=sqrt(1/beta(2)^2-pupmult3contour(i)^2);
     eta3(i)=sqrt(1/beta(3)^2-pupmult3contour(i)^2);
     eta4(i)=sqrt(1/beta(4)^2-pupmult3contour(i)^2);
     transmission(i)=(2*(mu(2)*eta2(i))/(mu(2)*eta2(i)+mu(1)*eta1(i)));
     reflection(i)=((mu(1)*eta1(i)-mu(2)*eta2(i))/(mu(2)*eta2(i)+mu(1)*eta1(i)))^3;
     PI4(i)=reflection(i)*transmission(i);
     dp(i)=(pupmult3contour(i)-pupmult3contour(i-1))/dt;
     Jupmult31(i)=imag(PI4(i)*pupmult3contour(i)^1.5*dp(i)/eta2(i));
     Jupmult32(i)=imag(PI4(i)*pupmult3contour(i)^0.5*dp(i));
end %for

for i=2:length(pdownref1contour)
    %downref1
     eta1(i)=sqrt(1/beta(1)^2-pdownref1contour(i)^2);
     eta2(i)=sqrt(1/beta(2)^2-pdownref1contour(i)^2);
     eta3(i)=sqrt(1/beta(3)^2-pdownref1contour(i)^2);
     eta4(i)=sqrt(1/beta(4)^2-pdownref1contour(i)^2);
     transmission(i)=(2*mu(2)*eta2(i))./(mu(1)*eta1(i)+mu(2)*eta2(i));
     reflection(i)=(mu(2)*eta2(i)-mu(3)*eta3(i))./(mu(2)*eta2(i)+mu(3)*eta3(i));
     if isreal(reflection(i))==0
         reflection(i)=real(reflection(i))+sqrt(-1)*-1*imag(reflection(i));
     end
     PI5(i)=transmission(i)*reflection(i);
     dp(i)=(pdownref1contour(i)-pdownref1contour(i-1))/dt;
     Jdownref11(i)=imag(PI5(i)*pdownref1contour(i)^1.5*dp(i)/eta2(i));
     Jdownref12(i)=imag(PI5(i)*pdownref1contour(i)^0.5*dp(i));
end %for

for i=2:length(pdownref2contour)
    %downref2
     eta1(i)=sqrt(1/beta(1)^2-pdownref2contour(i)^2);
     eta2(i)=sqrt(1/beta(2)^2-pdownref2contour(i)^2);
     eta3(i)=sqrt(1/beta(3)^2-pdownref2contour(i)^2);
     eta4(i)=sqrt(1/beta(4)^2-pdownref2contour(i)^2);
     transmission(i)=(2*mu(2)*eta2(i))/(mu(3)*eta3(i)+mu(2)*eta2(i));
     transmission2(i)=(2*mu(3)*eta3(i))/(mu(2)*eta2(i)+mu(3)*eta3(i));
     transmission3(i)=(2*mu(2)*eta2(i))/(mu(1)*eta1(i)+mu(2)*eta2(i));
     reflection(i)=real((mu(3)*eta3(i)-mu(4)*eta4(i))/(mu(3)*eta3(i)+mu(4)*eta4(i)));   
     if isreal(reflection(i))==0
         reflection(i)=real(reflection(i))+sqrt(-1)*-1*imag(reflection(i));
     end
     PI6(i)=transmission(i)*reflection(i)*transmission2(i)*transmission3(i);
     dp(i)=(pdownref2contour(i)-pdownref2contour(i-1))/dt;
     Jdownref21(i)=imag(PI6(i)*pdownref2contour(i)^1.5*dp(i)/eta2(i));
     Jdownref22(i)=imag(PI6(i)*pdownref2contour(i)^0.5*dp(i));
end %for

for i=2:length(pupdown1contour)
    %updown1
     eta1(i)=sqrt(1/beta(1)^2-pupdown1contour(i)^2);
     eta2(i)=sqrt(1/beta(2)^2-pupdown1contour(i)^2);
     eta3(i)=sqrt(1/beta(3)^2-pupdown1contour(i)^2);
     eta4(i)=sqrt(1/beta(4)^2-pupdown1contour(i)^2);
     transmission(i)=(2*mu(2)*eta2(i))/(mu(1)*eta1(i)+mu(2)*eta2(i));
     reflection1(i)=real((mu(2)*eta2(i)-mu(1)*eta1(i))/(mu(2)*eta2(i)+mu(1)*eta1(i)));   
     reflection2(i)=real((mu(2)*eta2(i)-mu(3)*eta3(i))/(mu(2)*eta2(i)+mu(3)*eta3(i)));   
     if isreal(reflection1(i))==0
         reflection1(i)=real(reflection1(i))+sqrt(-1)*-1*imag(reflection1(i));
     end
     if isreal(reflection2(i))==0
         reflection2(i)=real(reflection2(i))+sqrt(-1)*-1*imag(reflection2(i));
     end
     PI7(i)=transmission(i)*reflection1(i)*reflection2(i);
     dp(i)=(pupdown1contour(i)-pupdown1contour(i-1))/dt;
     Jupdown11(i)=imag(PI7(i)*pupdown1contour(i)^1.5*dp(i)/eta2(i));
     Jupdown12(i)=imag(PI7(i)*pupdown1contour(i)^0.5*dp(i));
end %for

for i=2:length(pupdown2contour)
    %updown2
     eta1(i)=sqrt(1/beta(1)^2-pupdown2contour(i)^2);
     eta2(i)=sqrt(1/beta(2)^2-pupdown2contour(i)^2);
     eta3(i)=sqrt(1/beta(3)^2-pupdown2contour(i)^2);
     eta4(i)=sqrt(1/beta(4)^2-pupdown2contour(i)^2);
     transmission1(i)=(2*mu(2)*eta2(i))/(mu(3)*eta3(i)+mu(2)*eta2(i));
     transmission2(i)=(2*mu(3)*eta3(i))/(mu(3)*eta3(i)+mu(2)*eta2(i));
     transmission3(i)=(2*mu(2)*eta2(i))/(mu(1)*eta1(i)+mu(2)*eta2(i));
     reflection1(i)=real((mu(2)*eta2(i)-mu(1)*eta1(i))/(mu(2)*eta2(i)+mu(1)*eta1(i)));   
     reflection2(i)=real((mu(3)*eta3(i)-mu(4)*eta4(i))/(mu(3)*eta3(i)+mu(4)*eta4(i)));   
     if isreal(reflection1(i))==0
         reflection1(i)=real(reflection1(i))+sqrt(-1)*-1*imag(reflection1(i));
     end
     if isreal(reflection2(i))==0
         reflection2(i)=real(reflection2(i))+sqrt(-1)*-1*imag(reflection2(i));
     end
     PI8(i)=transmission1(i)*transmission2(i)*transmission3(i)*reflection1(i)*reflection2(i);
     dp(i)=(pupdown2contour(i)-pupdown2contour(i-1))/dt;
     Jupdown21(i)=imag(PI8(i)*pupdown2contour(i)^1.5*dp(i)/eta2(i));
     Jupdown22(i)=imag(PI8(i)*pupdown2contour(i)^0.5*dp(i));
end %for

% figure(3)
% subplot(2,5,1), plot(PI1)
% subplot(2,5,2), plot(PI2)
% subplot(2,5,3), plot(PI3)
% subplot(2,5,4), plot(PI4)
% subplot(2,5,5), plot(PI5)
% subplot(2,5,6), plot(PI6)
% subplot(2,5,7), plot(PI7)
% subplot(2,5,8), plot(PI8)
% 
% figure(4)
% subplot(8,1,1), plot(xx,Jdirect1)
% subplot(8,1,2), plot(xx,Jupmult11)
% subplot(8,1,3), plot(xx,Jupmult21)
% subplot(8,1,4), plot(xx,Jupmult31)
% subplot(8,1,5), plot(xx,Jdownref11)
% subplot(8,1,6), plot(xx,Jdownref21)
% subplot(8,1,7), plot(xx,Jupdown11)
% subplot(8,1,8), plot(xx,Jupdown21)

J1=Jdirect1+Jupmult11+Jupmult21+Jupmult31+Jdownref11+Jdownref21+Jupdown11+Jupdown21;
J2=Jdirect2+Jupmult12+Jupmult22+Jupmult32+Jdownref12+Jdownref22+Jupdown12+Jupdown22;

con1=conv(rt,J1).*dt;
con2=conv(rt,J2).*dt;
dcon1=diff(con1)./dt;
dcon2=diff(con2)./dt;
sourcecon1=conv(dcon1,sthnorm).*dt;
sourcecon2=conv(dcon2,sthnorm).*dt;
v4(1,2)=0;
v4(1,1)=0;
v5(1,2)=0;
v5(1,1)=0;
v(1,2)=0;
v(1,1)=0;

for i=2:length(xx)
    v4(i,2)=A4*scale*sourcecon1(i);
    v4(i,1)=xx(i);
end
for i=2:length(xx)
    v5(i,2)=A5*scale*sourcecon2(i);
    v5(i,1)=xx(i);
end
for i=2:length(xx)
    v(i,1)=xx(i);
    v(i,2)=(v4(i,2)+v5(i,2));
end

%output v
return