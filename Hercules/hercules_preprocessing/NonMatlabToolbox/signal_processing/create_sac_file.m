%
% create sac file 
%
%   y=create_sac_file(start_time,dt,y,name)
%

function y=create_sac_file(start_time,dt,y,name)

N=length(y);
header_balin=newSacHeader(N,dt,start_time);
sacfile_balin=[name '.sac'];
wtSac(sacfile_balin,header_balin,y);

end