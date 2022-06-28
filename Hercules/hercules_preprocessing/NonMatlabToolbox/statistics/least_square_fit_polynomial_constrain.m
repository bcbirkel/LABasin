function [m] = least_square_fit_polynomial_constrain(dataX,dataY,polynomialconstrain)
%
% least_square_fit_polynomial_constrain: 
%
%   computes the least square fit asuming the abscense of any of the 
%   elements of a polynomial of order=length(polynomialconstrain)-1 where
%   polynomialconstrain is an array of 0 and 1 s that indicates if you have
%   an element or not, e.g.
%   
%   Linear fit
%    polynomialconstrain = [1 1]
%   Linear fit with constant value = 0
%    polynomialconstarin = [1 0]
%                                          

% Now it is just working with Linear fit with constant value = 0

m=sum(dataY)/sum(dataX);

y=m*dataX;
plot(dataX,y)

return
