%
% [stfinal]=computebandfilter_sac(signal,dt,fupper,poles,flower,passes)
%

function [stfinal]=computebandfilter_sac(signal,dt,fupper,poles,flower,passes,taperoffon)
     machine='n';
    n=length(signal);
    signal0=[signal];% 0*signal];
    
    create_sac_file(0,dt,signal0,'velocity');
    fp=fopen('process','w');        
    fprintf(fp,'%s\n',['read ' 'velocity.sac']);
    fprintf(fp,'%s\n','rmean');
    fprintf(fp,'%s\n','rtrend');   
    if(taperoffon==1)
        fprintf(fp,'%s\n','taper');
    end       
%    fprintf(fp,'%s\n',    ['hp n ' num2str(poles) ' corner ' num2str(flower) ' p '  num2str(passes)])  ;
%    fprintf(fp,'%s\n',    ['lp n ' num2str(poles) ' corner ' num2str(fupper) ' p '  num2str(passes)])  ; 
    fprintf(fp,'%s\n',    ['bp n ' num2str(poles) ' corner ' num2str(flower) ' ' num2str(fupper) ' p '  num2str(passes)])  ;
    fprintf(fp,'%s\n','write over')  ;
    
   
    fprintf(fp,'%s\n','quit');
    fclose(fp); 
    system( 'chmod 777 process');
    system( '/usr/local/sac/bin/sac < process ');
            
    % now read the signal again
    [t,st,SAChdr]= fget_sac('velocity.sac',0,machine);
    stfinal0=st.signal;
    stfinal0=st.signal;%fixbaseline(stfinal0,3);
    stfinal=stfinal0(1:n);
    system( 'rm -f process');
return

