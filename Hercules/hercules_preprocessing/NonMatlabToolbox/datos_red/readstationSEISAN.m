%**************************************************************************
%
%  Read a station using the SEISAN FORMAT
%
%  Author: Leonardo Ramirez-Guzman
%          II, UNAM, Mexico
%
%          Part of this function is based on Ana Laura Ruiz' (also at the
%          II UNAM) c function and the obspy.seisan.core python  code from
%          the ObsPy Python toolbox (Beyreuther et al, 2010, SRL electronic
%          seismologist)
%
%
%**************************************************************************
%clear
%close all
%pathName='/home/LRamirezG/Documents/March20EarthquakeMexico/processseisan/seisan_files/'
%fname=[pathName '2012-03-20-1802-43S.Regio_132'];
% %seisan=readstationSEISAN(fname)
 %
%earthquakes_to_review
%for iEarth=1:10
%fname=earthquake(iEarth).path
%fname='/home/LRamirezG/Documents/2012-06-29-1300-20S.MANII_066'%earthquake(1).path
function seisan=readstationSEISAN(fname);


stationsDatosRed=loaddatosred();

% Check whether a file is SIESAN or not and the byte order, architecture
% and version
fp=fopen(fname,'rb');
data=fread(fp,12*80,'char=>char');
seisan.byteorder_arch_version=seisan_getVersion(data);
fclose(fp);


% First line
fp=fopen(fname,'rb');
dataLine=seisan_readline(fp,80)';

seisan.networkname=dataLine(2:30);
seisan.numberofchannels=str2num(dataLine(31:33));
seisan.year  =str2num(dataLine(34:36))+1900;
seisan.doy   =str2num(dataLine(38:40));
seisan.month =str2num(dataLine(42:43));
seisan.day   =str2num(dataLine(45:46));
seisan.hr    =str2num(dataLine(48:49));
seisan.min   =str2num(dataLine(51:52));
seisan.sec   =str2num(dataLine(54:59));
seisan.totaltimewindow=str2num(dataLine(61:69));
                                    % 70-80 FREE

numberOfLines=  seisan.numberofchannels/3 + mod(seisan.numberofchannels, 3);
numberOfLines1=round(numberOfLines);
cort=0;
corre=0;
if (numberOfLines-numberOfLines1>.01)
    numberOfLines=numberOfLines1-1;
    cort=1;
    corre=round((numberOfLines1-numberOfLines)*3);
end
                                    
% Second line
dataLine=seisan_readline(fp,80)';

% Third line to the number of channel. In this implementation it is assumed
% that 3 channels are assigned per station
iChn=1;
%seisan=[]
iI=[2 28 54 ];
fI=[9 35 61];

for iSt=1:numberOfLines
   
    dataLine=seisan_readline(fp,80)';  
    
    for iChn=1:3
        seisan.data.station(iSt).channel(iChn).name=dataLine(iI(iChn):fI(iChn));
        seisan.data.station(iSt).component(iChn).starttimereleventfile=...
            str2num(dataLine(11:17));
        seisan.data.station(iSt).component(iChn).dataintervallength   =...
            str2num(dataLine(19:26));    
    end
    
    if(strcmp(seisan.data.station(iSt).channel(1).name(1:4),seisan.data.station(iSt).channel(2).name(1:4)));
        if(strcmp(seisan.data.station(iSt).channel(2).name(1:4),seisan.data.station(iSt).channel(3).name(1:4)));
            seisan.data.station(iSt).name=seisan.data.station(iSt).channel(1).name(1:4);           
            iChn=3+iChn;
        end
        
    else
       error('Stations do not have 3 consecutive channels')
    end            
end
if cort==1
    dataLine=seisan_readline(fp,80)';  
end
if (numberOfLines < 10)
        for i=1:10-numberOfLines
         seisan_readline(fp,80)';
        end
end

% Check that all stations are in DatosRed
for iSt=1:numberOfLines

    for iStDR=1:length(stationsDatosRed);
        if(strmatch(seisan.data.station(iSt).name,stationsDatosRed(iStDR).name(1:4),'exact'));
            stationsDatosRed(iStDR).name;
            seisan.data.station(iSt).lon=stationsDatosRed(iStDR).lon;
            seisan.data.station(iSt).lat=stationsDatosRed(iStDR).lat;
            seisan.data.station(iSt).ele=stationsDatosRed(iStDR).ele;
            seisan.data.station(iSt).factor=stationsDatosRed(iStDR).factor;
        end
end
end


%figureR
% Start reading channels and data
iSt=1;
for iChn=1:seisan.numberofchannels-corre
      
    dataLine=seisan_readline(fp,1040)';
    iComp=mod(iChn,3);
    if(iComp==0)
        iComp=3;
    end
    %seisan.station(iChn).component(A5=bString(:)
    seisan.data.station(iSt).component(iComp).filechannelname=[dataLine(1:4) dataLine(6:9)];

    seisan.data.station(iSt).component(iComp).century =str2num(dataLine(10)); 
    if(seisan.data.station(iSt).component(iComp).century == 0);
        century=1900;
    else
        century=2000; 
    end
    
    
    seisan.data.station(iSt).component(iComp).year = century+str2num(dataLine(11:12));
    seisan.data.station(iSt).component(iComp).doy  = str2num(dataLine(14:16));
    seisan.data.station(iSt).component(iComp).month= str2num(dataLine(18:19));
    seisan.data.station(iSt).component(iComp).day  = str2num(dataLine(21:22));
    seisan.data.station(iSt).component(iComp).hr   = str2num(dataLine(24:25));
    seisan.data.station(iSt).component(iComp).min  = str2num(dataLine(27:28));
    seisan.data.station(iSt).component(iComp).sec  = str2num(dataLine(30:35));

    seisan.data.station(iSt).component(iComp).samplingrate  = str2num(dataLine(37:43));
    seisan.data.station(iSt).component(iComp).numsamples    = str2num(dataLine(44:50));
    seisan.data.station(iSt).component(iComp).gainfactorS   = dataLine(76);
    if(strcmp(seisan.data.station(iSt).component(iComp).gainfactorS,'G')==0)
       seisan.data.station(iSt).component(iComp).gainfactor.gainfactor=1;
    else
        seisan.data.station(iSt).component(iComp).gainfactor = str2num(dataLine(148:159));
    end
    
    
   seisan.data.station(iSt).component(iComp).intsize       = str2num(dataLine(77));
    if(isempty(seisan.data.station(iSt).component(iComp).intsize));
         seisan.data.station(iSt).component(iComp).intsize    =2     ;
    end

    seisan.data.station(iSt).component(iComp).signal=[];

    if( seisan.data.station(iSt).component(iComp).intsize == 4);
        stdata=fread(fp,seisan.data.station(iSt).component(iComp).numsamples+2,'int32');
        seisan.data.station(iSt).component(iComp).signal=stdata(3:length(stdata));
    end
    if( seisan.data.station(iSt).component(iComp).intsize == 2);
        stdata=fread(fp,seisan.data.station(iSt).component(iComp).numsamples+2,'int16');   
        seisan.data.station(iSt).component(iComp).signal=stdata(3:length(stdata)-10);
    end
    if(iComp==3);
       iSt=iSt+1;
    end  
     %figureR;
     %signal=seisan.data.station(iSt).component(iComp).signal
     %numSamples=length(signal)
     %meanRm=mean(signal(4:round(numSamples/10)));
     %plot(signal(3:length(signal)-20));
%     title( seisan.station(iChn).filechannelname);
     %drawnow

end
fclose(fp);






%end