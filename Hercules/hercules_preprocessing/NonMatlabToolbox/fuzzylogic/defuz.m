% a=defuz(x,y) DEFUZZIFICA LA FUNCION Y EN EL DOMINIO X
% UTILIZANDO EL DEFIZZIFICADOR POR CENTROIDE DE ÁREA: EL 
% CENTROIDE ES a.
function c=defuz(y,mu)

mul=y.*mu;
num=sum(mul);
den=sum(mu);

c=num/den;