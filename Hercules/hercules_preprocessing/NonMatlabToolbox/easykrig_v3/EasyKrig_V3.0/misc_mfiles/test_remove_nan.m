function		[x,y,z,v]=test_remove_nan(filename)
%%% remove NaN's
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

dat=load(filename);
x=dat(:,2);
y=dat(:,1);
if size(dat,2) >= 4   % 3D data
   z=dat(:,3);
   v=dat(:,4);
   indx=find( isnan(x) | isnan(y) | isnan(z) | isnan(v));
   x(indx)=[];
   y(indx)=[];
   z(indx)=[];
   v(indx)=[];
else                    % 2D data
   v=dat(:,3);
   indx=find( isnan(x) | isnan(y) | isnan(v));
   x(indx)=[];
   y(indx)=[];
   v(indx)=[];
   z=[];
end