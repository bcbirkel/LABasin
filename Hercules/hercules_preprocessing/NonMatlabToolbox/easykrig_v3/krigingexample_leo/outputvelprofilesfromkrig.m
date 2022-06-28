
% Script to output a mat file with the important info from the krigging
% processing


vel=data.out.krig.Vg
nx=data.out.krig.nx
ny=data.out.krig.ny
nz=data.out.krig.nz

xx=data.out.krig.Xg;
yy=data.out.krig.Yg;
zz=data.out.krig.Zg;

save examplevel.mat vel nx ny nz xx yy zz

