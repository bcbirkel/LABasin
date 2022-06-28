%**************************************************************************
%
%  Read a station using the BMSF
%
%**************************************************************************

%clear
%pathName='/home/lramirezguzman/Desktop/RSM/data/11_12_2011/'
%fname=[pathName 'CCT31112.111'];


function station = readstationBMSF(fname)

[header, data] = hdrload(fname);
% Go through the header to get the location of the station
sizeHeader = size(header);

% *******LOCALIZACION DE LA ESTACION **************************************

stringLocationSt='CLAVE DE LA ESTACION';
iCount=0;
doit=1;
while (doit==1)
    iCount=iCount+1;
    if (isempty(strfind(header(iCount,:), stringLocationSt)))
        doit=1;
    else
        doit=0;
        % The String was found thus retrieve the requiered info
        strfind(header(iCount,:), stringLocationSt)
        [token , remain ] = strtok(header(iCount,:),':');
        [token1, remain1] = strtok(remain);
        [token2, remain2] = strtok(remain1);
        station.code = token2;
    end
end

stringLocationSt='COORDENADAS DE LA ESTACION';
iCount=0;
doit=1;
while (doit==1)
    iCount=iCount+1;
    if (isempty(strfind(header(iCount,:), stringLocationSt)))
        doit=1;
    else
        doit=0;
        % The String was found thus retrieve the requiered info
        strfind(header(iCount,:), stringLocationSt)
        [token , remain ] = strtok(header(iCount,:),':');
        [token1, remain1] = strtok(remain);
        [token2, remain2] = strtok(remain1);
        station.lat=str2num(token2);
        clear('token','token1','token2','remain','remain1','remain2')
        
        strfind(header(iCount+1,:), stringLocationSt)
       [token , remain ] = strtok(header(iCount+1,:),':');
        [token1, remain1] = strtok(remain);
        [token2, remain2] = strtok(remain1);
        station.lon=str2num(token2);
        clear('token','token1','token2','remain','remain1','remain2')
        strfind(header(iCount+1,:), stringLocationSt)
       [token , remain ] = strtok(header(iCount+2,:),':');
        [token1, remain1] = strtok(remain);
        [token2, remain2] = strtok(remain1);
        station.elev=str2num(token2);
    end
end

stringLocationSt='COORDENADAS DEL EPICENTRO';
iCount=0;
doit=1;
while (doit==1)
    iCount=iCount+1;
    if (isempty(strfind(header(iCount,:), stringLocationSt)))
        doit=1;
    else
        doit=0;
        % The String was found thus retrieve the requiered info
        strfind(header(iCount,:), stringLocationSt)
        [token , remain ] = strtok(header(iCount,:),':');
        [token1, remain1] = strtok(remain);
        [token2, remain2] = strtok(remain1);
        station.source_lat=str2num(token2);
        clear('token','token1','token2','remain','remain1','remain2')
        strfind(header(iCount+1,:), stringLocationSt)
        [token , remain ] = strtok(header(iCount+1,:),':');
        [token1, remain1] = strtok(remain);
        [token2, remain2] = strtok(remain1);
        station.source_lon=str2num(token2);
        clear('token','token1','token2','remain','remain1','remain2')
        strfind(header(iCount+1,:), stringLocationSt)
        [token , remain ] = strtok(header(iCount+2,:),':');
        [token1, remain1] = strtok(remain);
        [token2, remain2] = strtok(remain1);
        station.source_depth=str2num(token2);
    end
end

stringLocationSt='HORA EPICENTRO (GMT)';
iCount=0;
doit=1;
while (doit==1)
    iCount=iCount+1;
    if(isempty(strfind(header(iCount,:), stringLocationSt)))
        doit=1;
    else
        doit=0;
        % The String was found thus retrieve the requiered info
        strfind(header(iCount,:), stringLocationSt)
        [token , remain ] = strtok(header(iCount,:),':');
        [token1, remain1] = strtok(remain);
        [token2, remain2] = strtok(remain1);
        station.source_hor=str2num(token2(1:2));
        station.source_min=str2num(token2(4:5));
        station.source_sec=str2num(token2(7:end));
    end
end


%******************CHANNEL ORIENTATION ************************************
clear stringToSearch
stringToSearch='NUMERO DE CANALES';
iCount=0;
doit=1;
while (doit==1)
    iCount=iCount+1;
    if (isempty(strfind(header(iCount,:), stringToSearch)))
        doit=1;
    else
        doit=0;
        % The String was found thus retrieve the requiered info
        strfind(header(iCount,:), stringToSearch)
        [token , remain ] = strtok(header(iCount,:),':');
        [token1, remain1] = strtok(remain);
        [token2, remain2] = strtok(remain1);
        station.numberOfChannels=str2num(token2);
        clear('token','token1','token2','remain','remain1','remain2')
    end
end


stringToSearch='CANAL-';
iCount=0;
doit=1;
while (doit==1)
    iCount=iCount+1;
    if (isempty(strfind(header(iCount,:), stringToSearch)))
        doit=1;
    else
        doit=0;
        % The String was found thus retrieve the requiered info
        strfind(header(iCount,:), stringToSearch)
        remain=header(iCount+1,:);
        for iChannel=1:station.numberOfChannels
            [token , remain1 ] = strtok(remain,' ');
            clear remain;
            remain=remain1;
            station.channel(iChannel).string=token;
            clear token;
        end
    end
end

% Now obtain the data
for iCh=1:station.numberOfChannels
    station.channel(iCh).acc = data(:,iCh);
end

%******************Delta t per channel ************************************
stringLocationSt='INTERVALO DE MUESTREO, C1-C6';
iCount=0;
doit=1;
while (doit==1)
    iCount=iCount+1;
    if (isempty(strfind(header(iCount,:), stringLocationSt)))
        doit=1;
    else
        doit=0;
        % The String was found thus retrieve the requiered info
        strfind(header(iCount,:), stringToSearch)
        remain=header(iCount,:);
        for iChannel=1:station.numberOfChannels
            [token , remain1 ] = strtok(remain,'/');
            clear remain;
            remain=remain1;
            [token , remain1 ] = strtok(remain,'/');
            if(isempty(token))
                station.channel(iChannel).dt=(station.channel(iChannel-1).dt);
            else
             station.channel(iChannel).dt = str2num(token);
            end
            clear token;            
        end                      
    end
end


stringLocationSt='HORA DE LA PRIMERA MUESTRA (GMT)';
iCount=0;
doit=1;
while (doit==1)
    iCount=iCount+1;
    if (isempty(strfind(header(iCount,:), stringLocationSt)))
        doit=1;
    else
        doit=0;
        % The String was found thus retrieve the requiered info
        strfind(header(iCount,:), stringLocationSt)
        [token , remain ] = strtok(header(iCount,:),':');
        [token1, remain1] = strtok(remain,':');
        station.origin_hor=str2num(token1);
        [token2, remain2] = strtok(remain1,':');
        station.origin_min=str2num(token2);
        [token3, remain3] = strtok(remain2,':');
        station.origin_sec=str2num(token3);
    end
end

% %% Display each channel
% figureR
% for iP=1:station.numberOfChannels
%     subplot(station.numberOfChannels,1,iP)
%     plot(station.channel(iP).acc)
%     grid
% end

return
end
