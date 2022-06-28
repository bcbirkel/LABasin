% Read Waveform data, correct instrument response, filter, decimate, and
% attached processed data to Eqk Database.
% You may want to check the files since PZ files from IRIS may fuck you up.

close all; clear; clc;

% Load event information
load EarthquakeDatabase.mat

% Frequencies for filtering, Time and Window
FMin = 1/100;
FMax = 1;
Dt = 0.25;
Time = (0:Dt:200).';
NT = length(Time);
W = tukeywin(NT,0.25);

% Save Parameters
save Parameters.mat FMin FMax Dt Time W;

% Go and find waveform files by identifying SACPZ
for nS = 1:size(Eqk,1)
    % Find all SACPZ of the SCEDC stations and read station parameters
    % In the previous step we potentially deleted what isn't useful
    PZ_Files = dir(['/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/ObsData/Velocity/SACPZ/SACPZ*']);
    
    % Allocate memory
    R_Ntw = cell(size(PZ_Files));
    R_Code = cell(size(PZ_Files));
    R_Loc = cell(size(PZ_Files));
    R_Chn = cell(size(PZ_Files));
    R_Dsc = cell(size(PZ_Files));
    R_Lat = cell(size(PZ_Files));
    R_Lon = cell(size(PZ_Files));
    R_Elev = cell(size(PZ_Files));
    R_Depth = cell(size(PZ_Files));
    R_Azm = cell(size(PZ_Files));
    R_Dip = cell(size(PZ_Files));
    
    % Open short.sitechan file
    Fid1 = fopen(['/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/hercules_preprocessing/MatlabCodes/AuxFiles/' ...
        num2str(Eqk(nS).Year) num2str(Eqk(nS).Month,'%02d') num2str(Eqk(nS).Day,'%02d') ...
        num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') ...
        '.sitechan'],'w');
    
    % Open StationChanel.txt
    Fid3 = fopen(['/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/hercules_preprocessing/MatlabCodes/AuxFiles/' ...
        num2str(Eqk(nS).Year) num2str(Eqk(nS).Month,'%02d') num2str(Eqk(nS).Day,'%02d') ...
        num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') ...
        '_STN.txt'],'w');
    
    % Save short.site
     Fid2 = fopen(['/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/hercules_preprocessing/MatlabCodes/AuxFiles/' ...
        num2str(Eqk(nS).Year) num2str(Eqk(nS).Month,'%02d') num2str(Eqk(nS).Day,'%02d') ...
        num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') ...
        '.site'],'w');
    
    for nf = 1:size(PZ_Files,1)
        % Open SACPZ File
        Fid = fopen(['/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/ObsData/Velocity/SACPZ/' PZ_Files(nf).name]);
        
        % Read SACPZ File line by line
        fgetl(Fid);         % First Line Just *s
        
        NL = fgetl(Fid);    % Network
        In = 23;
        if length(NL) < 23
            In = 17;
        end
        R_Ntw{nf} = NL(In:end);
        
        NL = fgetl(Fid);    % Station
        R_Code{nf} = NL(In:end);
        
        NL = fgetl(Fid);    % Location 00, 10,...
        R_Loc{nf} = NL(In:end); R_Loc{nf} = strrep(R_Loc{nf},'  ','--');
        
        NL = fgetl(Fid);    % Channel
        R_Chn{nf} = NL(In:end);
        
        fgetl(Fid);         % CREATED
        
        fgetl(Fid);         % START
        
        fgetl(Fid);         % END
        
        NL = fgetl(Fid);    % Description
        R_Dsc{nf} = NL(In:min(In+49,length(NL)));
        
        NL = fgetl(Fid);    % Latitude
        R_Lat{nf} = NL(In:end);
        
        if isnan(str2double(R_Lat{nf}))
            error(['No Lat ' PZ_Files(nf).name])
        end
        NL = fgetl(Fid);    % Longitude
        R_Lon{nf} = NL(In:end);
        
        if isnan(str2double(R_Lon{nf}))
            error(['No Lon ' PZ_Files(nf).name])
        end
        
        NL = fgetl(Fid);    % Elevation
        R_Elev{nf} = NL(In:end);
        
        if isnan(str2double(R_Elev{nf}))
            error(['No Elev ' PZ_Files(nf).name])
        end
        
        NL = fgetl(Fid);	% Depth
        R_Depth{nf} = NL(In:end);
        
        if isnan(str2double(R_Depth{nf}))
            error(['No Depth ' PZ_Files(nf).name])
        end
        
        NL = fgetl(Fid);    % Dip
        R_Dip{nf} = NL(In:end);
        
        NL = fgetl(Fid);    % Azimuth
        R_Azm{nf} = NL(In:end);
        

        fclose(Fid);

        fprintf(Fid1,'%-6s %-9s 2000001       -1       -1 n     %8.4f %6.1f %6.1f - %+57s 00:00:00 \n',...
            R_Code{nf},R_Chn{nf},str2double(R_Depth{nf}),str2double(R_Azm{nf}),str2double(R_Dip{nf}),'01/01/00');
        
        fprintf(Fid3,'%-6s %-8s %-7s %8.4f %10.4f %10.4f %10.4f %6s %8.1f %8.1f\n',...
            R_Ntw{nf},R_Code{nf},R_Chn{nf},str2double(R_Lat{nf}),str2double(R_Lon{nf}),...
            str2double(R_Elev{nf}),str2double(R_Depth{nf}),R_Loc{nf},str2double(R_Azm{nf}),str2double(R_Dip{nf}));
    end
    
    % Find unique STN
    [~,In,~] = unique(R_Code,'stable');
    In(end+1) = size(R_Code,1);
    
    % Look at all the stations
    for nr = 1:size(In,1)-1
        
        % Obs will have this size
        Wave = zeros(In(nr+1)-In(nr),NT);
                
        % Match SACPZ and SAC waveform file
        for nc = 1:In(nr+1)-In(nr)
            if R_Ntw{In(nr)+nc-1} == 'CI'
                SAC_Path = ['/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/ObsData/Velocity/' ...
                    R_Ntw{In(nr)+nc-1} '.' R_Code{In(nr)+nc-1} '.' R_Loc{In(nr)+nc-1} '.' R_Chn{In(nr)+nc-1} ...
                    '.D.' num2str(Eqk(nS).Year)  '.' num2str(Eqk(nS).JDay,'%03d') '.' ...
                    num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') '.SAC'];
                SAC_File = dir(SAC_Path);
            
            elseif R_Ntw{In(nr)+nc-1} == 'CE'
                SAC_Path = ['/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/ObsData/Velocity/' ...
                    R_Ntw{In(nr)+nc-1} '_' R_Code{In(nr)+nc-1} '_BH'  R_Chn{In(nr)+nc-1}(3) '_vel.SAC'];
                SAC_File = dir(SAC_Path);
                
            else 
                SAC_Path = '';
            end
            
            
            if isfile(SAC_Path)
               
                % Load observed seismogram
               Sac = rsac('little-endian',SAC_Path);
               
               if size(Sac,1) == NT
                   
                % Data time Sampling and window
                DtOld = Sac(1,3);
                
                % Decimate
                Sac = decimate(rtrend(Sac(:,2)),round(Dt/DtOld));
                
                % Correct instrument response
                Wave(nc,:) = transfer(W.*rtrend(Sac(1:NT)),NT,1/Dt,['/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/currentSim/ObsData/Velocity/SACPZ/' ...
                        PZ_Files(In(nr)+nc-1).name],FMin,FMax);
%                else
%                    fprintf('bad file rip')
               end
            end
        end
                    
                    
         % Create Structure for station
        Eqk(nS).Stn(nr).Network = R_Ntw{In(nr)};
        Eqk(nS).Stn(nr).Name = R_Code{In(nr)};
        Eqk(nS).Stn(nr).Level = R_Loc{In(nr)};
        Eqk(nS).Stn(nr).Description = R_Dsc{In(nr)};
        Eqk(nS).Stn(nr).Lat = str2double(R_Lat{In(nr)});
        Eqk(nS).Stn(nr).Lon = str2double(R_Lon{In(nr)});
        Eqk(nS).Stn(nr).Elev = str2double(R_Elev{In(nr)});
        Eqk(nS).Stn(nr).Depth = str2double(R_Depth{In(nr)});
    
        % Create at least three components
        Channels = {'BHZ','BHN','BHE'};
        Azims = [0 0 90];
        Dips = [-90 0 0];
        for nc1 = 1:3*ceil((In(nr+1)-In(nr))/3)
            Eqk(nS).Stn(nr).Channels(nc1).Name =  Channels{1+mod(nc1,3)};
            Eqk(nS).Stn(nr).Channels(nc1).Data = zeros(1,NT);
            Eqk(nS).Stn(nr).Channels(nc1).Azimuth = Azims(1+mod(nc1,3));
            Eqk(nS).Stn(nr).Channels(nc1).Dip = Dips(1+mod(nc1,3));
        end
        
        for nc1 = 1:In(nr+1)-In(nr)
            if  strcmp(R_Chn{In(nr)+nc1-1},'BHN')
                Eqk(nS).Stn(nr).Channels(1).Name =  R_Chn{In(nr)+nc1-1};
                Eqk(nS).Stn(nr).Channels(1).Data = Wave(nc1,:);
                Eqk(nS).Stn(nr).Channels(1).Azimuth = str2double(R_Azm{In(nr)+nc1-1});
                Eqk(nS).Stn(nr).Channels(1).Dip = str2double(R_Dip{In(nr)+nc1-1});
            elseif strcmp(R_Chn{In(nr)+nc1-1},'BHE')
                Eqk(nS).Stn(nr).Channels(2).Name =  R_Chn{In(nr)+nc1-1};
                Eqk(nS).Stn(nr).Channels(2).Data = Wave(nc1,:);
                Eqk(nS).Stn(nr).Channels(2).Azimuth = str2double(R_Azm{In(nr)+nc1-1});
                Eqk(nS).Stn(nr).Channels(2).Dip = str2double(R_Dip{In(nr)+nc1-1});
            elseif strcmp(R_Chn{In(nr)+nc1-1},'BHZ')
                Eqk(nS).Stn(nr).Channels(3).Name =  R_Chn{In(nr)+nc1-1};
                Eqk(nS).Stn(nr).Channels(3).Data = Wave(nc1,:);
                Eqk(nS).Stn(nr).Channels(3).Azimuth = str2double(R_Azm{In(nr)+nc1-1});
                Eqk(nS).Stn(nr).Channels(3).Dip = str2double(R_Dip{In(nr)+nc1-1});
            end
        end
        
        fprintf(Fid2,'%-7s 2000001       -1  %8.4f %9.4f  %8.4f %-50s -    -         0.0000    0.0000 01/01/00 00:00:00 \n',...
            R_Code{In(nr)},str2double(R_Lat{In(nr)}),str2double(R_Lon{In(nr)}),0.001*str2double(R_Elev{In(nr)}),R_Dsc{In(nr)});
    end
    fclose all;
end

fprintf('A4 finished \n')

% Save Earthquake Structure
save EarthquakeDatabase.mat Eqk;
