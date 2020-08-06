clear all; clc; format long g;
%-------------------------------------------------------------------------%
% CU-PSHA Software                                                        %
%                                                                         %
%                                                                         %
%   Evaluate percent of occurrence of ground shaking (MMI) of interest    %
%                                                                         %
%                                                                         %
%                                                                         %
%-------------------------------------------------------------------------%

%(Step 1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load neccesary data(Require) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

MMI            = 7;   % Percent of exceedance
Yr             = 50;  % Year              
MG             = 13;  
% Remarks for MG is to determine MMI vs PGA(unit g) relationship 
% 1. Cancani (1904)	                    Log(PGA) = 0.33MMI-1.17
% 2. Ishimoto (1932)	                Log(PGA) = 0.4228MMI-1.4178
% 3. Kawasumi (1951)	                Log(PGA) = 0.3409MMI-0.5705
% 4. Neumann (1954)	                    Log(PGA) = 0.3041MMI-0.0096
% 5. Savarensky and Kirnos (1955)	    Log(PGA) = 0.2773MMI-0.5915
% 6. Hershberger (1956)	                Log(PGA) = 0.4329MMI-0.921
% 7. Richter (1958)	                    Log(PGA) = 0.3297MMI-0.499
% 8. Medvedev and Sponheuer (1969)	    Log(PGA) = 0.3019MMI-0.2394
% 9. Okamoto (1973)  	                Log(PGA) = 0.32MMI-0.3985
% 10. Trifunac and Brady (1975)-Vert.	Log(PGA) = 0.304MMI-0.2115
% 11. Trifunac and Brady (1975)-Horiz.	Log(PGA) = 0.2949MMI+0.046
% 12. Theodulides and Papazachos (1992)	Log(PGA) = 0.291MMI+0.304
% 13. Shabestari and Yamazaki (2001)	Log(PGA) = 0.2545MMI+0.2977

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert from m/sec2 to g
if      MG     == 1
PGA            = (10^((0.33*MMI)-1.17))/981;       % 1. Cancani (1904)
else if MG     == 2
PGA            = (10^((0.4228*MMI)-1.4178))/981;   % 2. Ishimoto (1932)
else if MG     == 3
PGA            = (10^((0.3409*MMI)-0.5705))/981;   % 3. Kawasumi (1951)
else if MG     == 4
PGA            = (10^((0.3041*MMI)-0.0096))/981;   % 4. Neumann (1954)
else if MG     == 5
PGA            = (10^((0.2773*MMI)-0.5915))/981;   % 5. Savarensky and Kirnos (1955)
else if MG     == 6
PGA            = (10^((0.4329*MMI)-0.921))/981;    % 6. Hershberger (1956)
else if MG     == 7
PGA            = (10^((0.3297*MMI)-0.499))/981;    % 7. Richter (1958)
else if MG     == 8
PGA            = (10^((0.3019*MMI)-0.2394))/981;   % 8. Medvedev and Sponheuer (1969)
else if MG     == 9
PGA            = (10^((0.32*MMI)-0.3985))/981;     % 9. Okamoto (1973)
else if MG     == 10
PGA            = (10^((0.304*MMI)-0.2115))/981;    % 10. Trifunac and Brady (1975)-Vert.
else if MG     == 11
PGA            = (10^((0.29493*MMI)+0.046))/981;   % 11. Trifunac and Brady (1975)-Horiz.
else if MG     == 12
PGA            = (10^((0.291*MMI)+0.304))/981;     % 12. Theodulides and Papazachos (1992)
else if MG     == 13
PGA            = (10^((0.2545*MMI)+0.2977))/981;   % 13. Shabestari and Yamazaki (2001)
end
end
end
end
end
end
end
end
end
end
end
end
end

%(Step 2)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Call the hazard-curve file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HZC1           = fopen('Output/Hazard-curve.txt','r'); 
[HZC2,HZCcount]= fscanf(HZC1,'%f',[202 inf]);
HZCpoint       = HZCcount/202;
HZC            = HZC2';
fclose(HZC1);
Acc            = 0.005:0.01:1.995; % in g unit (1g = 9.81 m/s2)

%(Step 3)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Print out the percent map file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ske1           = strcat('Output/MMI',num2str(MMI),'Y', num2str(Yr),'.txt');
Ske2           = fopen(Ske1,'wt');
fprintf(Ske2,'%s\t%s\t%s\t\n' ,['longitude ','latitude ','Probability(%)']);

%(Step 4)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculate the percent map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rate1          = zeros(HZCpoint-1,200);
for i          = 1:HZCpoint-1
for j          = 1:200
Rate1(i,j)     = HZC(i+1,j+2); 
end
end

for i          = 1:HZCpoint-1
Rate           = Rate1(i,:);
for k          = 2:200

if PGA<Acc(k)  & PGA>Acc(k-1)       
MM(i)          = interp1(Acc(k-1:k),Rate(k-1:k),PGA,'spline');
else if PGA    ==Acc(k)  
MM(i)          = Rate(k);
else if PGA    ==Acc(1)
MM(i)          = Rate(1);
else if PGA    > Acc(200)
MM(i)          = 0;
else if  PGA   < Acc(1)
MM(i)          = 1;
end    
end
end
end
end
end

MMGG(i) =(1-exp((-1)*MM(i)*Yr))*100;
%MMGG(i) =1/MM(i);
fprintf(Ske2,'\n%f %f %f',[HZC(i+1,1);HZC(i+1,2);MMGG(i)]);

end
fclose(Ske2);

%(Step 5)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Show status (finish) of percent map calculation %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

disp([' '])
disp(['********************************************************************************* '])
disp(['  Finish calculation of percent of groundshaking (MMI) occurence'])
disp(['  Percent of ground shaking equal or greater than MMI ' num2str(MMI) ' in ' num2str(Yr) ' year'])
disp(['********************************************************************************* '])
disp([' '])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%