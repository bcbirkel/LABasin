function	initialization()
% function initialization perform initialization

global para color data hdl

    
curDir=pwd;							% current directory
para.curDir=curDir;
cmd=['which(''startkrig'')'];
startfilename=eval(cmd);
indx=find(startfilename == '.');
HDIR=startfilename(1:max(indx)-11);
DirBackslash=length(find(HDIR == '\'));
DirColon=length(find(HDIR == ':'));
DirSlash=length(find(HDIR == '/'));
HDIRpath=['' HDIR ''];

if DirBackslash >= 1			% WINDOWS
   DirMark='\';
   AddPathCmd=['addpath '  '''' HDIRpath '\general''', ...
                ' ''' HDIRpath '\dataprep''', ' ''' HDIRpath '\variogram''', ...
                ' ''' HDIRpath '\krig''', ' ''' HDIRpath '\visualization''', ...
                ' ''' HDIRpath '\help''', ' ''' HDIRpath '\images''', ...
                ' ''' HDIRpath '\bin''', ' ''' HDIRpath '\data''', ...
                ' ''' HDIRpath '\output''', ' ''' HDIRpath '\''',  ' -begin'];
   PlatForm=1;
   markersize=6;
elseif DirSlash >= 1			% UNIX/LINUX
   DirMark='/';
   AddPathCmd=['addpath '  '''' HDIRpath '/general''', ...
                ' ''' HDIRpath '/dataprep''', ' ''' HDIRpath '/variogram''', ...
                ' ''' HDIRpath '/krig''', ' ''' HDIRpath '/visualization''', ...
                ' ''' HDIRpath '/help''', ' ''' HDIRpath '/images''', ...
                ' ''' HDIRpath '/bin''', ' ''' HDIRpath '/data''', ...
                ' ''' HDIRpath '/output''', ' ''' HDIRpath '/''',  ' -begin'];
   PlatForm=2;  
   markersize=8;
end
eval(AddPathCmd)

%% Optimization Toolbox
OptimOption=0;
pp=path;
pindx=find(pp == DirMark);
if pindx+5 <= length(pp)
   ppl=length(pindx);
else
   ppl=length(pindx)-1;
end
for i=1:ppl
   StrIndx=pindx(i)+1:pindx(i)+5;
   if max(StrIndx) <= length(pp)
     DirStr=pp(StrIndx);
     if DirStr == 'optim'
       OptimOption=1;
     end
   end
end

data.in.dim=2;              % default is 2D case

% set parameters
para.home_dir=HDIR;
para.platform=PlatForm;
para.optim=OptimOption;
para.status=0;
para.file_dir.data_conversion=HDIR;
para.file_dir.datafile=HDIR;
para.file_dir.data_format_file=HDIR;
para.file_dir.gridfile=HDIR;
para.file_dir.parafile=HDIR;
para.file_dir.batch_filename=HDIR;
para.file_dir.batch_log=HDIR;
para.file_dir.mat_file_in=HDIR;
para.file_dir.mat_file_out=HDIR;

para.dataprep.filename='';
para.dataprep.ext_prog=0;
para.dataprep.dat_conv_fname='';
para.dataprep.xy_switch=0;
para.krig.data_format_file=[];
para.status.dataprepfig=0;
para.status.dataprep=0;
para.status.variogramfig=0;
para.status.variogram=0;
para.status.krigingfig=0;
para.status.kriging=0;
para.status.dispkrigfig=0;
para.status.dispkrig=0;

hdl.status.dataprepfig=0;
hdl.status.variogramfig=0;
hdl.status.krigingfig=0;
hdl.status.dispkrigfig=0;

para.dataprep.ytox=1;
para.dataprep.ztox=1;
para.dataprep.ext_prog=0;
para.dataprep.filter_type=2;				% default filter = mean
para.dataprep.reduct_fac=1;
para.dataprep.filter_supt=1;
para.krig.load_data_format_file=0;
para.dataprep.data_disptype=1;
para.dataprep.data_disptype0=1;

para.vario.max_nugt=1;
para.vario.max_sill=1.5;
para.vario.max_powr=4.0;
para.vario.max_range=sqrt(2);			% normalized range
para.vario.max_hole=4*pi/para.vario.max_range;
para.vario.max_lscl=para.vario.max_range;
para.vario.load_para=0;
para.vario.para_file='';
para.krig.load_para=0;
para.krig.vario_para=0;
para.krig.krig_para=0;
para.krig.both_para=1;
para.krig.para_file_in='';
para.krig.load_data_file=0;
para.krig.batch_file_proc=0;
para.krig.bat_proc_cnt=0;
para.krig.data_file='';

para.dispkrig.markersize=markersize;
para.dispkrig.customized_grid_data_markersize=4;
para.krig.load_griddata_file=0;
para.dispkrig.trackline.dispflag=1;

% color
color.background=[0.8 0.8 0.8];
color.grey=[0.75 0.75 0.75];
color.dark_grey=[0.65 0.65 0.65];
color.blue=[0.753 0.753 0.753];

% data structure
data.in.dim=2;          % default

% graphic handle
hdl.object.edit_w=0.04;
hdl.object.pushbtn_w=0.05;
hdl.object.pushbtn_l=0.14;
hdl.object.popmenu=0.04;
hdl.object.radio_w=0.03;
hdl.object.txt_w8=0.03;
hdl.object.txt_w10=0.05;

if para.platform == 2   % If Unix OS
   getXvisual;
end
hdl.msg.h0=[];

warning off

