% Load DatosRed.txt
%fname='DatosRed.txt'
%[header, data]=hdrload(fname);
function stationsDatosRed=loaddatosred()
fp=fopen('DatosRed.txt','r');
numStations= str2num(fgetl(fp));
for i=1:numStations
   lineTmp=fgetl(fp);
   [token , remain ] = strtok(lineTmp,',');
   stationsDatosRed(i).name  =token;
   clear token;
   [token , remain ] = strtok(remain,',');   
   stationsDatosRed(i).lon   =str2num(token);
   clear token;
   [token , remain ] = strtok(remain,',');   
   stationsDatosRed(i).lat   =str2num(token);
   clear token;
   [token , remain ] = strtok(remain,',');   
   stationsDatosRed(i).ele   =str2num(token);
   clear token;
   [token , remain ] = strtok(remain,',');   
   stationsDatosRed(i).factor   =str2num(token);
end
