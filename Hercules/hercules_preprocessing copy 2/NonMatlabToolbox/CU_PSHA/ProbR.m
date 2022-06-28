function PR    = ProbR(long,lat,Fnum)
format long g;
%-------------------------------------------------------------------------%
% CU-PSHA Software                                                        %
%                                                                         %
%                                                                         %
%    Calculate probability of distance from earthquake source to site     %
%                                                                         %
%                                                                         %
%                                                                         %
%-------------------------------------------------------------------------%

%(Step 1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Call EQ source %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qn1            = strcat('Input/',num2str(Fnum),'.txt');
Qn2            = fopen(Qn1,'r');
[Qn3,Qncount]  = fscanf(Qn2,'%f',[3 inf]); 
Qnpoint        = Qncount/3;
Qn             = sortrows(Qn3');
fclose(Qn2);

%(Step 2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Determine probability of distance distribution %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rdist          = zeros(Qnpoint);
for i          = 1:Qnpoint
Rdist(i)       = sqrt((((Qn(i,1)-long)*100)^2)+(((Qn(i,2)-lat)*100)^2)+(Qn(i,3)^2)); % distance is in km unit
end
distmax        = max(Rdist);
distmin        = min(Rdist);
  
subR           = distmin+((distmax-distmin)/100):(distmax-distmin)/50:distmax-((distmax-distmin)/100);
DistR          = hist(Rdist,subR);
probdistR      = DistR/sum(DistR);

ProbR          = zeros(2,50);
ProbR(1,:)     = subR;
ProbR(2,:)     = probdistR; 
PR             = ProbR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remarks:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Supplementary scripts for print out the probability of distance %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Distname     = strcat('Output/Prob-distance (Long ',num2str(long),' Lat ',num2str(lat),' Fault no. ',num2str(Fnum),').txt');
% PD           = fopen(Distname,'wt');
% fprintf(PD,'%s\t%s\t' ,['Distance',' Probability', ]);
% fprintf(PD,'\n%f %f',[PR]);
% fclose(PD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
