clear

%--------------------------------------------------------------------------
% smooth krig velocity
%--------------------------------------------------------------------------

load examplevel.mat

% Put a limit to the lowest value (you decide)
k=find(vel<3);
vel(k)=3;

% smoothing (read smooth3 for more info)
velS=smooth3(vel,'gaussian',[11 11 11],4);

figureR; 
slice(xx,yy,-zz,velS,[-92 -90 -88 ],37,[-10])
xlabel('Lon','FontSize',14)
ylabel('Lat','FontSize',14)
shading interp
loni=xx(1,:,1);
lati=yy(:,1,1);
z=reshape(zz(1,1,:),nz,1);

%--------------------------------------------------------------------------
% NOW OBTAIN THE PROFILES A
%--------------------------------------------------------------------------
%PROFILES
fp=fopen('example.fun','w')
numProfiles=length(loni)*length(lati);
[lo,la]=meshgrid(loni,lati);
longitude=reshape(lo,length(loni)*length(lati),1);
latitude =reshape(la,length(loni)*length(lati),1);

% NUMBER OF VELOCITY PROFLES
figureR;
hold
fprintf(fp,'\n %d ',numProfiles)
ilonlat=0;
for iLat=1:length(lati)
    for iLon=1:length(loni)
        ilonlat=ilonlat+1;
        
        for i=1:length(velS(1,1,:));
            %CHECK THAT YOU DO NOT HAVE NAN, IF YOU DO REPEAT THE PREVIOUS
            %VALUE
            if( isnan(velS(iLat,iLon,i))==0 )
                vs(i)=velS(iLat,iLon,i) ;
            else
                if(i == 1)
                    vs(i)=min(velS(iLat,iLon,:));
                    if(isnan(vs(i))==1)
                        vs(i)=3.5;
                    end                   
                    
                else
                    vs(i)=vs(i-1);
                end
            end
            
            [vpBr,vsBr,rhoBr]=computebrochersproperties(0,vs(i),6.14);
            vs(i)=vsBr/1000;
            vp(i)=vpBr;
            rho(i)=rhoBr;
        end
        %       zInt = linspace(prunc(iLat,iLon),60,40);
        zInt =  linspace(0,200,60);
        vsInt =  interp1(z,vs,zInt);
        vpInt =  interp1(z,vp,zInt);
        rhoInt = interp1(z,rho,zInt);
        
        
        
        %IF NOT DEFINED THE VELOCITY REPEAT IT
        for i=1:length(vsInt)
            if(i>1)
                if( isnan(vsInt(i))==1 )
                    vsInt(i)=vsInt(i-1) ;
                    vpInt(i)=vpInt(i-1) ;
                    rhoInt(i)=rhoInt(i-1) ;
                end
            else
                if(i==1)
                    if( isnan(vsInt(i))==1 )
                        vsInt(i)=min(vsInt(:)) ;
                        vpInt(i)=min(vpInt(:)) ;
                        rhoInt(i)=min(rhoInt(:)) ;
                    end
                end
            end
        end
               
        
        % WRITE NUMBER OF POINTS IN THE PROFILE, ONE DEPTH AND PROFILE PER
        % LINE
        %        zRef=(zInt-prunc(iLat,iLon))*1000;
        zRef=(zInt);
        fprintf(fp,'\n %d',length(zRef))
        fprintf(fp,'\n ')
        fprintf(fp,' %e',zRef);
        fprintf(fp,'\n ')
        fprintf(fp,' %e',vpInt); %Vp
        fprintf(fp,'\n ')
        fprintf(fp,' %e',vsInt*1000); %Vs
        fprintf(fp,'\n ')
        fprintf(fp,' %e',rhoInt); %Density
        plot(vsInt,zRef);axis ij
        clear vp;
        clear vs;
        clear rho;
    end
end


fp=fopen('example.surf','w')
count=0;
fprintf(fp,' %d %d \n ', length(lati), length(loni));
fprintf(fp,' %e ',lati)

fprintf(fp,' \n ');
fprintf(fp,' %e ',loni)
fprintf(fp,' \n ');

for iLon=1:length(loni)
    fprintf(fp,'\n');
    for iLat=1:length(lati)
        fprintf(fp,' %d ', count);
        count=count+1;
    end
    
end
fclose(fp);

