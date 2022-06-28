function plotdata(x,data,varargin)

%PLOTDATA  displays seismograms as a horizontal record section
%
%  plotdata(data) displays seismograms given by the columns of data ordered 
%  from first to last on the y-axis centered on the corresponding event number.
%
%  plotdata(data,scale) normalizes each record by the scale factor (scalar).
%
%  plotdata(data,'overlay') adds an additional plot after the last event of 
%  all the traces superposed, which is useful for comparing repeating events.
%
%  plotdata(data,scale,'overlay')  performs all of the above.

%  David Schaff  8-8-98
%  could modify to allow scale to be a vector and also a 'clipping' flag

switch nargin 
 case 2
   scale = 1;
   OVERLAY = 0;
   line = 'b';
 case 3
   if ischar(varargin{1})
      if strcmp('overlay', varargin{1})
         OVERLAY = 1;
         line = 'b';
      else
         OVERLAY = 0;
         line = varargin{1};
      end
      scale = 1;
   else
      OVERLAY = 0;
      scale = varargin{1};
      line = 'b';
   end
 case 4
   if ischar(varargin{1})
      if strcmp('overlay', varargin{1})
         OVERLAY = 1;
         if ischar(varargin{2})
            line = varargin{2};
            scale = 1;
         else
            line = 'b';
            scale = varargin{2};
         end
      else
         line = varargin{1};
         if ischar(varargin{2}) & strcmp('overlay', varargin{2})
            OVERLAY = 1;
            scale = 1;
         else
            OVERLAY = 0;
            scale = varargin{2};
         end
      end
   else
      scale = varargin{1};
      if ischar(varargin{2})
         if strcmp('overlay', varargin{2}) 
            OVERLAY = 1;
            line = 'b';
         else
            OVERLAY = 0;
            line = varargin{2};
         end
      end
   end
 case 5
   OVERLAY = 1;
   if ischar(varargin{1})
      if strcmp('overlay', varargin{1})
         if ischar(varargin{2})
            line = varargin{2};
            scale = varargin{3};
         else
            line = varargin{3};
            scale = varargin{2};
         end
      else
         line = varargin{1};
         if ischar(varargin{2}) & strcmp('overlay', varargin{2})
            scale = varargin{3};
         else
            scale = varargin{2};
         end
      end
   else
      scale = varargin{1};
      if ischar(varargin{2})& strcmp('overlay', varargin{2}) 
         line = varargin{3};
      else
         line = varargin{2};
      end
   end
end
   
% useful for clipped data
data=data*diag(1./mean(abs(data)));
stack=mean(data,2);
stackmax=max(stack);
data=.65*data/stackmax;

data=data/scale;
data=-data;		% change polarity for plotting b/c y-axis reversed

[numSamp,numEve]=size(data);
separate=ones(numSamp,1)*(1:numEve);		
plot(x, data+separate,line)
addone=1;

%plot(t(:,ones(numEve,1)),data+shift,'b') could do later for t vector.

if OVERLAY
   hold on
   together=(numEve+1)*ones(size(data));		
%   plot(data+together,'b')
   plot(x, data+together)
   hold off
   addone = 1 + addone;
end

% make plot pretty
ylim=get(gca,'ylim');
axis tight
set(gca,'ydir','reverse','ylim',[0 numEve+addone])
box on
