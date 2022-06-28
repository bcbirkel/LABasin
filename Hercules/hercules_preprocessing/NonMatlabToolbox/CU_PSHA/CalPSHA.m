function Pacc  = CalPSHA(long,lat)
format long g;
%-------------------------------------------------------------------------%
% CU-PSHA Software                                                        %
%                                                                         %
%                                                                         %
%                      Calculate hazard curve                             %
%                                                                         %
%                                                                         %
%                                                                         %
%-------------------------------------------------------------------------%

%(Step 1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load earthquake source parameters(Require) %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
QP1              = fopen('Input/EQ-parameters.txt','r');
[QP2,QPcount]    = fscanf(QP1,'%f',[10 inf]);
Fnum             = QPcount/10;
QP               = sortrows(QP2');
fclose(QP1);
HZZ              = zeros(Fnum,200); %

%(Step 2)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Determine Probability of magnitude and distance %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for f            = 1:Fnum
PM               = ProbM(QP(f,1),QP(f,3),QP(f,4),QP(f,5),QP(f,6),QP(f,7),QP(f,8),QP(f,9)); %(10 case study)
PR               = ProbR(long,lat,QP(f,1)); %(50 case study)

%(Step 3)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Determine mean and standard deviation of Peaks Ground Accerelation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lnPHA            = zeros(10,50); 
SDPHA            = zeros(10,50); 
for i            = 1:10
for j            = 1:50
if QP(f,10)      == 1
PGA              = att1(PM(i,1),PR(1,j),QP(f,2));
else if QP(f,10) == 2
PGA              = att2(PM(i,1),PR(1,j),QP(f,2));
else if QP(f,10) == 3
PGA              = att3(PM(i,1),PR(1,j),QP(f,2));
else if QP(f,10) == 4
PGA              = att4(PM(i,1),PR(1,j),QP(f,2));
else if QP(f,10) == 5
PGA              = att5(PM(i,1),PR(1,j),QP(f,2));
else if QP(f,10) == 6
PGA              = att6(PM(i,1),PR(1,j),QP(f,2));
else if QP(f,10) == 7
PGA              = att7(PM(i,1),PR(1,j),QP(f,2));
end
end
end
end
end
end
end
lnPHA(i,j)       = PGA(1); % mean PGA in unit "ln(gal)" (1 gal=0.01 m/s2)
SDPHA(i,j)       = PGA(2); % standard deviation
end
end

%(Step 4)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Evaluate probability of exceedance in each PGA level, M, and R %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remarks : 1) Set maximum ground shaking level of interest = 2g
%           2) Set PGA interval range = 0.01g (200 case study from 0-2g)  

PHA              = 0.005:0.01:1.995; % in g unit (1g = 9.81 m/s2)
for k            = 1:200 
lnacc(k)         = log(PHA(k)*981);  % "Change g unit to be ln(gal)unit"
for i            = 1:10
for j            = 1:50
probPHA(k,i,j)   = 1-normcdf(lnacc(k),lnPHA(i,j),SDPHA(i,j));
end
end
end

%(Step 5)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Evaluate the average rate of threshold magnitude exceedance %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delM2            = 0.5;              
MC               = QP(f,3)-delM2;
for k            = 1:200
for i            = 1:10
for j            = 1:50
if PM(i,1)       >= MC
lamda(k,i,j)     = PM(11,1)*PM(i,2)*PR(2,j)*probPHA(k,i,j);   
else 
lamda(k,i,j)     = PM(11,2)*PM(i,2)*PR(2,j)*probPHA(k,i,j); 
end
end
lamdasum1(k,i)   = sum(lamda(k,i,:));
end
lamdasum2(k)     = sum(lamdasum1(k,:)); %(200x1)
end
HZZ(f,:)         = lamdasum2(1,:)';
end
HZS              = sum (HZZ,1);

%(Step 6)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Sum the probability of exceedance from every EQ sources %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pat              = zeros(1,202);
Pat(1,1)         = long;
Pat(1,2)         = lat;
for i            = 1:200
Pat(1,i+2)       = HZS(i); %(1x2002)
end
Pacc             = Pat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
