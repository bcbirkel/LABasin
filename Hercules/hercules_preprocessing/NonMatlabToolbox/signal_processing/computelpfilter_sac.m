%
% [stfinal]=computebandfilter_sac(signal,dt,fupper,poles,flower,passes)
%

function [stfinal]=computelpfilter_sac(signal,dt,fupper,poles,passes)
    machine='n';

    create_sac_file(0,dt,signal,'velocity');
    fp=fopen('process','w');
    
    
    fprintf(fp,'%s\n',['read ' 'velocity.sac']);
    fprintf(fp,'%s\n','taper');
    fprintf(fp,'%s\n','rmean');
    fprintf(fp,'%s\n','rtrend');
    fprintf(fp,'%s\n',['lp n ' num2str(poles) ' corner ' num2str(fupper) ' p '  num2str(passes)])  ;
    fprintf(fp,'%s\n','taper');
    fprintf(fp,'%s\n','write over')  ;
    fprintf(fp,'%s\n','quit');
    fclose(fp); 
    system( 'chmod 777 process');
    system( 'sac < process ');
            
    % now read the signal again
    [t,st,SAChdr]= fget_sac('velocity.sac',0,machine);
    stfinal=st.signal;
    system( 'rm -f process');
return

