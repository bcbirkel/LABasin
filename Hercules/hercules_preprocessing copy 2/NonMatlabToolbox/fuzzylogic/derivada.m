
function [df]=derivada(dx,fx)
%Nombre dl programa:"Diferencias Finitas"
% [df]=derivada(dx,fx) calcula la derivada de la funcion fx evaluda en
% puntos con una separacion dx usando la fórmula de diferencias adelantadas 
% f: Valores de la funcion que se desea derivar
% dx: diferencial de x
%
% Variables de salida:
% df: vector con las derivadas en cada punto

M = length(fx); 

% Derivada con diferencias finitas centradas
for i=1:(M-1)
    dfxc(i)=(fx(i+1)-fx(i))/(dx);
end
% En los extremos se calcula con diferencias hacia adelante y hacia atras;
    dfxc(M)=(fx(M)-fx(M-1))/dx;
df=dfxc';
    

