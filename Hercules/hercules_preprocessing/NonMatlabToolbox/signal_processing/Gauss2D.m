%EFILTER constructs a band pass Gaussian filter
%   E = EFILTER(Size, Sigma1, Sigma2) designs a 2-dimensional Gaussian
%   digital filter where Size is a two element 
%   [rows, cols] vector specifying the size of the filter to construct, 
%   Sigma1 and Sigma2, the normalized standard deviation on the major 
%   and minor axes. E is normalized to have max(E)=1;
%
%   If E = EFilter(Size, Sigma1, Sigma2, alpha), where alpha is an angle
%   in radians, it will return and plot an Gaussian filter rotated
%   counter-clockwise through alpha.
%
%   If E = EFilter(Size, Sigma1, Sigma2, alpha, xoff, yoff), where xoff
%   and yoff are offsets (mu1 and mu2) in the x and y direction, 
%   it will return and
%   plot an Gaussian filter which is offset by the specified amount.
%   An offset of 0 corresponds to the center and an offset of 1
%   corresponds to the edge of the filter. A positive offset shifts the
%   filter in the positive direction.
%
%   Calling EFilter(...) without assigning the output variable
%   plots the 3D surface described by the function.

% Katie Streit   kstreit@rice.edu
% ELEC 301
% Rice University
% Alan Juarez alanjuar@usc.edu
% ZHS-266
% University of Southern California
%
% December 2001 - Jan 2017

% Much of this code was based on Peter Kovesi's  (pk@cs.uwa.edu.au)
% Matlab function for a lowpass Butterworth filter.


function varargout = Gauss2D(sze, S1, S2, varargin);

if nargin == 3
alpha = 0;
offx = 0;
offy = 0;
elseif nargin == 4
offx = 0;
offy = 0;
alpha = varargin{1};
elseif nargin == 6
alpha = varargin{1};
offx = varargin{2};
offy = varargin{3};
else
	error('Invalid number of input arguments');
end

    if nargout > 1
        error('Invalid number of output arguments');
    end
    

%%extracts the sizes from sze 
rows = 2*sze(1)-1;
cols = 2*sze(2)-1;

%x and y matrices normalized to +/-.5 and an offset of offx or offy
[x,y]=meshgrid(([1:cols]-fix(cols/2))/cols-offx/2,...
    ([1:rows]-fix(rows/2))/rows-offy/2);

%applies a linear transformation to rotate through alpha. Note that it takes
% uses negative alpha, which is caused by x and y being independent matrices.
x2 = (x*cos(alpha) - y*sin(-alpha));
y2 = (x*sin(-alpha) + y*cos(alpha));

%constructs Gaussian 
exponent = ((x2/S1).^2 + (y2/S2).^2);
f       = (exp(-exponent));  


f = f/max(max(f));

m=ceil(rows/2);
p=ceil(cols/2);

if nargout > 0
  varargout{1} = f(m:end,p:end);
else
  %Plots a normalized (+/- 1), interpolated 3D image of the filter
  surf([-1:2/(cols-1):1],[-1:2/(rows-1):1], f);
  shading interp;
  title('Elliptical Butterworth filter');
  xlabel('x');
  ylabel('y');
  zlabel('intensity');
  grid on;
end
