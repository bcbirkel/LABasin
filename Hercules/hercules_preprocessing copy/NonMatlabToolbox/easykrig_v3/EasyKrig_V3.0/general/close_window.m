function			close_window(window_index)
% funstion			close_window(window_index)
% close window based on the window index
%%
%%  Kriging Software Package  version 3.0,   May 1, 2004
%%  Copyright (c) 1999, 2001, 2004, property of Dezhang Chu and Woods Hole Oceanographic
%%  Institution.  All Rights Reserved.

global hdl 

switch window_index
	case 1			% data preparation window
      close(hdl.dataprep.h0);
      hdl.status.dataprepfig=0;
	case 2			% semi-variogram/correlogram window
      close(hdl.vario.h0);
      hdl.status.variogramfig=0;
   case 3			% kriging	window
      close(hdl.krig.h0);
      hdl.status.krigingfig=0;
   case 4		% visulization window
      close(hdl.dispkrig3d.h0);
      hdl.status.dispkrigfig=0;
      if hdl.status.krigingfig == 1
         p=findobj(hdl.krig.h0,'type','axes');
         if ~isempty(p) delete(p);end
         return
      end
   case 5		% 2D-3D variogram/correlogram visulization window
      close(hdl.dispvario2d3d.h0);
      hdl.status.dispvario2d3d=0;
end
   
   
return