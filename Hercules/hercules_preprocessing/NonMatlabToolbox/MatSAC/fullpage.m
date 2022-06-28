function props=fullpage(off)

%FULLPAGE  resizes figure on screen and paper to fit entire page.
%
%  FULLPAGE resizes the current figure so that it will print on
%     the entire page of an 8.5 x 11 piece of paper and resizes
%     the window to appear how it will print.
%  FULLPAGE OFF returns these figure properties to their default
%     settings.
%  FULLPAGE ON is equivalent to FULLPAGE.
%
%  properties=FULLPAGE; returns a cell with the fullpage properties
%  that can be used as figure(properties{:}).  Current figure is not
%  modified however.
%
%  See also FIGURE, ORIENT, FIGTITLE, FIGCAPTION

%  David Schaff 8-21-98

if nargout==1
  props={'PaperPosition',[0.25    0.6    8   9.65], ...
         'Units','inches','Position',[0 -1    8   9.65]};

elseif nargin == 0 | strcmp(off,'on')
   set(gcf,'PaperPosition',[0.25    0.6    8   9.65])
   set(gcf,'Units','inches','Position',[0   -1    8   9.65])
elseif strcmp(off,'off')
   set(gcf,'Units',get(0,'DefaultFigureUnits'), ...
	'Position',get(0,'DefaultFigurePosition'), ...
	'PaperPosition',get(0,'DefaultFigurePaperPosition'))
else
   disp('To return to default settings type:  fullpage off.')

end	% end nargout
