function   yout=bessel0(x)
%% Cylindrical Bessel function of order 1

x0=5;
[m,n]=size(x);
indx1=find(x < x0);
indx2=find(x >= x0);
x1=x(indx1)/3;
x2=x(indx2);
c1=2.2499997;c2=1.2656208;c3=.3163866;c4=0.0444479;c5=0.0039444;
c6=0.00021;
y1=1-c1*x1.^2+c2*x1.^4-c3*x1.^6+c4*x1.^8-c5*x1.^10+c6*x1.^12;
y2=sqrt(2./(pi*x2)).*cos(x2-pi/4);

y(indx1)=y1;
y(indx2)=y2;
yout=reshape(y,m,n);


%z=besselj(0,x);

%z-yout
%plot(x,z,x,y,'.-')
