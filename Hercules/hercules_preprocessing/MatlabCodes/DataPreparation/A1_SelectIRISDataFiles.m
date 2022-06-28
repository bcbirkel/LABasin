% Select files with good quality data:
% Estimate signal-to-noise ratio larger than 8

close all; clear; clc;

% Load event information
load EarthquakeDatabase.mat;

% Go and find waveform files by identifying SACPZ
for nS = 1:size(Eqk,1)
    
    % Find all SACPZ of the SCEDC stations and read station parameters
    PZ_Files = dir(['/home/scec-00/birkel/EarthquakePhysics/ObservedData/' ...
        num2str(Eqk(nS).Year) num2str(Eqk(nS).Month,'%02d') num2str(Eqk(nS).Day,'%02d') ...
        num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') ...
        '/SCEDC/SACPZ*']);
    
    % Allocate memory for PZ metadata
    R_Ntw = cell(size(PZ_Files));
    R_Code = cell(size(PZ_Files));
    R_Loc = cell(size(PZ_Files));
    R_Chn = cell(size(PZ_Files));
    R_Lat = cell(size(PZ_Files));
    R_Lon = cell(size(PZ_Files));
    
    % Loop over all the PZ files
    for nf = 1:size(PZ_Files,1)
        
        % Open SACPZ File
        Fid = fopen(['/home/scec-00/birkel/EarthquakePhysics/ObservedData/' ...
            num2str(Eqk(nS).Year) num2str(Eqk(nS).Month,'%02d') num2str(Eqk(nS).Day,'%02d') ...
            num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') ...
            '/SCEDC/' PZ_Files(nf).name]);
        
        % Read SACPZ Files line by line
        fgetl(Fid);         % First Line Just *'s
        
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
        fgetl(Fid);    % Description
        
        
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
        
        fgetl(Fid);    % Elevation
        fgetl(Fid);	% Depth
        fgetl(Fid);    % Dip
        fgetl(Fid);    % Azimuth
        
        fclose(Fid);
    end
    
    % Find unique STN
    [~,In,~] = unique(R_Code,'stable');
    In(end+1) = size(R_Code,2);
    
    % Look at all the stations
    for nr = 1:size(In,1)-1
        
        % Calculate station epicentral distance for P-Wave travel time
        [Dist,~,~] = vincentyinv(Eqk(nS).Lat,Eqk(nS).Lon,str2double(R_Lat{In(nr)}),str2double(R_Lon{In(nr)}));
        
        % Approximate P-wave traveltime from the prem model
        TauP = tauptime('mod','prem','dep',Eqk(nS).Depth,'km',Dist);
        TauP = TauP(1).time;
        
        % Match SACPZ and SAC waveform file
        for nc = 1:In(nr+1)-In(nr)
            
            % Find SAC waveform file
            SAC_File = dir(['/home/scec-00/birkel/EarthquakePhysics/ObservedData/' ...
                num2str(Eqk(nS).Year) num2str(Eqk(nS).Month,'%02d') num2str(Eqk(nS).Day,'%02d') ...
                num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') ...
                '/' R_Ntw{In(nr)+nc-1} '.' R_Code{In(nr)+nc-1} '.' R_Loc{In(nr)+nc-1} '.' R_Chn{In(nr)+nc-1} ...
                '.D.' num2str(Eqk(nS).Year)  '.' num2str(Eqk(nS).JDay,'%03d') '.' ...
                num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') '.SAC']);
            
            % If the SAC file doesn't exist, delete the PZ file
            if size(SAC_File,1) == 0
                disp(['File not found ' R_Ntw{In(nr)+nc-1} '.' R_Code{In(nr)+nc-1} '..' R_Chn{In(nr)+nc-1} ...
                    '.D.' num2str(Eqk(nS).Year)  '.' num2str(Eqk(nS).JDay,'%03d') '.' num2str(Eqk(nS).JDay,'%03d') '.' ...
                    num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') '.SAC']);
                
                system(['rm -r /home/scec-00/birkel/EarthquakePhysics/ObservedData/' ...
                    num2str(Eqk(nS).Year) num2str(Eqk(nS).Month,'%02d') num2str(Eqk(nS).Day,'%02d') ...
                    num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') ...
                    '/SCEDC/' PZ_Files(In(nr)+nc-1).name]);
                
            else
                % If file exists, load observed seismogram
                Sac = rsac('little-endian',['/home/scec-00/birkel/EarthquakePhysics/ObservedData/' ...
                    num2str(Eqk(nS).Year) num2str(Eqk(nS).Month,'%02d') num2str(Eqk(nS).Day,'%02d') ...
                    num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') ...
                    '/' R_Ntw{In(nr)+nc-1} '.' R_Code{In(nr)+nc-1} '.' R_Loc{In(nr)+nc-1} '.' R_Chn{In(nr)+nc-1} ...
                    '.D.' num2str(Eqk(nS).Year)  '.' num2str(Eqk(nS).JDay,'%03d') '.' ...
                    num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') '.SAC']);
                
                % Time Sampling
                Dt = Sac(1,3);
                
                % Waveform
                Sac = rtrend(Sac(:,2));
                
                % P-wave time sample
                P = round(TauP/Dt);
                
                % Window normalized Noise and Signal
                Noise = max(abs(Sac(1:P)));
                Signal = max(abs(Sac(P:end-1)));
                
                % Estimate signal to noise. If SNR<8 delete SAC and PZ files
                if Signal/Noise < 8
                    display(['Signal-to-noise of '  R_Ntw{In(nr)+nc-1} '.' R_Code{In(nr)+nc-1} '.' R_Loc{In(nr)+nc-1} '.' ...
                        R_Chn{In(nr)+nc-1} '.D.' num2str(Eqk(nS).Year)  '.' num2str(Eqk(nS).JDay,'%03d') '.' ...
                        num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') '.SAC' ' = ' ...
                        num2str(Signal/Noise)]);
                    
                    system(['rm -r /home/scec-00/birkel/EarthquakePhysics/ObservedData/' ...
                        num2str(Eqk(nS).Year) num2str(Eqk(nS).Month,'%02d') num2str(Eqk(nS).Day,'%02d') ...
                        num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') ...
                        '/' R_Ntw{In(nr)+nc-1} '.' R_Code{In(nr)+nc-1} '.' R_Loc{In(nr)+nc-1} '.' R_Chn{In(nr)+nc-1} ...
                        '.D.' num2str(Eqk(nS).Year)  '.' num2str(Eqk(nS).JDay,'%03d') '.' ...
                        num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') '.SAC']);
                    
                    system(['rm -r /home/scec-00/birkel/EarthquakePhysics/ObservedData/' ...
                        num2str(Eqk(nS).Year) num2str(Eqk(nS).Month,'%02d') num2str(Eqk(nS).Day,'%02d') ...
                        num2str(Eqk(nS).Hour,'%02d') num2str(Eqk(nS).Min,'%02d') num2str(Eqk(nS).Sec,'%02d') ...
                        '/SCEDC/' PZ_Files(In(nr)+nc-1).name]);
                end
            end
        end
    end
end