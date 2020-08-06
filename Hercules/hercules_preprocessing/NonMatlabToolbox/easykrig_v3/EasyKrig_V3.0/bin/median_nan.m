function    pmedian=median_nan(A,k)
%% function    pmean=median_nan(A,k)
%% computes the median value which ignores all nan's 
%% if A is an 1D array, k is not necessary, if A is a matrix
%% k is optional. Without k or k = 1, A is averaged over column,
%% and k = 2, average is over rows
%%
%%  Kriging Software Package  version 2.0,   October 29, 1999
%%  Copyright (c) 1999, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

D=size(A);
if D(1) == 1 | D(2) == 1						% 1-D array
   [indx]=find(~isnan(A));
   if ~isempty(indx)
      pmedian=median(A(indx));
   else
      pmedian=nan;
   end
else
   if nargin == 1
     k=1;				% default direction: average over each colume
  end
  if k == 1
    for i=1:D(2)
       [indx]=find(~isnan(A(:,i)));
       if ~isempty(indx)
          pmedian(i)=median(A(indx,i));
       else
          pmedian(i)=nan;
       end
    end
  else
    for i=1:D(1)
       [indx]=find(~isnan(A(i,:)));
       if ~isempty(indx)
          pmedian(i)=median(A(i,indx));
       else
          pmedian(i)=nan;
       end
   end
   pmedian=pmedian(:);
  end
end
  
