clear all; clc; format long g;
%-------------------------------------------------------------------------%
% CU-PSHA Software                                                        %
%                                                                         %
%                                                                         %
%            Evaluate the ground shaking of exceedance                    %
%                                                                         %
%                                                                         %
%                                                                         %
%-------------------------------------------------------------------------%

%(Step 1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load neccesary data(Require) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

Pc             = 2; % Percent of exceedance
Yr             = 50; % Year              
Lamda          = log10((-log(1-(Pc/100)))/Yr);

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
%%%% Print out the ground shaking map file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ske1           = strcat('Output/P',num2str(Pc),'Y',num2str(Yr),'.txt');
Ske2           = fopen(Ske1,'wt');
fprintf(Ske2,'%s\t%s\t%s\t\n' ,['longitude ','latitude ','PGA(g)']);

%(Step 4)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculate the ground shaking map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rate           = zeros(HZCpoint-1,200);
for i          = 1:HZCpoint-1
for j          = 1:200
if HZC(i+1,j+2)> 0
Rate(i,j)      = log10(HZC(i+1,j+2)); 
% take log10 into Annualsum because the values are extremely varies more than normal scale.
else
Rate(i,j)      = -19.3920890976347;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k          = 2:200
if Lamda>Rate(i,k) & Lamda<Rate(i,k-1)
PGA(i)         = interp1(Rate(i,k-1:k),Acc(k-1:k),Lamda,'spline');
else if Lamda  == Rate(i,k)
PGA(i)         = Acc(k);
else if Lamda  == Rate(i,1)
PGA(i)         = 0.005;
else if Lamda  > Rate(i,1)
PGA(i)         = 0;
else if Lamda  < Rate(i,200)
PGA(i)         = 2.000;
end
end
end
end
end
end
fprintf(Ske2,'\n%f %f %f',[HZC(i+1,1);HZC(i+1,2);PGA(i)]); 
end
fclose(Ske2);

%(Step 5)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Show status (finish) of ground shaking map calculation %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([' '])
disp(['********************************************************************************* '])
disp(['       Finish calculation of ground shaking'])
disp(['       ' num2str(Pc) ' % probability of exceedance in ' num2str(Yr) ' year'])
disp(['********************************************************************************* '])
disp([' '])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

