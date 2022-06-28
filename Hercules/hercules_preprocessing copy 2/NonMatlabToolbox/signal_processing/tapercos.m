function w = tapercos(v,l1,l2,f)
% w = tapercos(v,l1,l2,f)
% Tapered cosine (Tukey Window) function

x = (v-l1)/(l2-l1);
r = f;%*l1/(l2-l1);
w = 0*x;

for i = 1:length(x)
    if x(i) >= 0 && x(i) < r/2
        w(i) = 0.5*(1+cos(2*pi*(x(i)-r/2)/r));
    elseif x(i) >= r/2 && x(i) <= 1-r/2
        w(i) = 1;
    elseif x(i) >= 1-r/2 && x(i) <= 1
        w(i) = 0.5*(1+cos(2*pi*(x(i)-1+r/2)/r));
    end
end