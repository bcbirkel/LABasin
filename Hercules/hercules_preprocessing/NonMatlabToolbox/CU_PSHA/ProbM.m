function PM    = ProbM(Fnum,Mmax,Mmin,a,b,SR,Af,FMD)
format long g;
%-------------------------------------------------------------------------%
% CU-PSHA Software                                                        %
%                                                                         %
%                                                                         %
%            Calculate probability of earthquake magnitude                %
%                                                                         %
%                                                                         %
%                                                                         %
%-------------------------------------------------------------------------%

%(Step 1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Determine the required earthquake parameters %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
c                = 1.5;              % Hanks and Kanamori(1979) >> log(Mo=cM+d)
d                = 16.1;             % Hanks and Kanamori(1979) >> log(Mo=cM+d)
Mshear           = 3*(10^11);        % Rigidity or shear modulus (dyne/cm2)
alpha            = 2.3026*a;         % G-R relationship log(N)=a-bM -->> ln(N)=alpha-(beta*M) 
beta             = 2.3026*b;         % G-R relationship log(N)=a-bM -->> ln(N)=alpha-(beta*M) 
delM1            = 1;                % Young and Coppersmith (1985)
delM2            = 0.5;              % Young and Coppersmith (1985)
MC               = Mmax-delM2;       % Characteristic EQ (Young and Coppersmith, 1985)
Momax            = 10^((c*Mmax)+d);  % Maximum seismic moment(M0) 
fm               = zeros(10,1);

%(Step 2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Determine Probability Density Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subM             = (Mmin+((Mmax-Mmin)/20):(Mmax-Mmin)/10:Mmax-((Mmax-Mmin)/20))';
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1. Exponential magnitude distribution model (Youngs and Coppersmith, 1985) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FMD           == 1 
for i            = 1:10
fm(i)            = (((beta*exp((-beta)*(subM(i)-0)))/(1-exp((-beta)*(Mmax-0)))))*((Mmax-0)/10); 
end
Mull1            = (Mshear*Af*SR*(c-b)*(1-exp((-beta)*(Mmax-Mmin))))/(b*Momax*exp((-beta)*(Mmax-Mmin)));
Mull2            = exp(alpha-(beta*Mmin));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.2. Characteristic earthquake model (Youngs and Coppersmith, 1985; Convertito et al., 2006)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else if FMD      == 2 
Cat              =((beta*exp((-1)*(Mmax-0-delM1-delM2)))/(1-exp((-beta)*(Mmax-0-delM2))))*delM2;
for i            = 1:10
if subM(i)       < Mmin
fm(i)            = 0;
else if subM(i)  >= Mmin & subM(i)<=MC
fm(i)            = (((beta*exp((-beta)*(subM(i)-0)))/(1-exp((-beta)*(Mmax-0-delM2))))*(1/(1+Cat)))*((Mmax-0)/10);       
else if subM(i)  >= MC & subM(i)<=Mmax
fm(i)            = (((beta*exp((-beta)*(Mmax-0-delM1-delM2)))/(1-exp((-beta)*(Mmax-0-delM2))))*(1/(1+Cat)))*((Mmax-0)/10);  
else if subM(i)  > Mmax   
fm(i)            = 0;
end
end
end
end
end
K                = ((b*10^((-1)*c*delM2))/(c-b))+((b*exp(beta*delM1)*(1-(10^((-1)*c*delM2))))/c);
alphaNC          = (Mshear*Af*SR*(1-exp((-beta)*(Mmax-Mmin-delM2))))/(K*Momax*exp((-beta)*(Mmax-Mmin-delM2)));
Mull1            = (alphaNC*(beta*delM2*exp((-beta)*(Mmax-Mmin-delM1-delM2))))/(1-exp((-beta)*(Mmax-Mmin-delM2))) ;
Mull2            = exp(alpha-(beta*Mmin));
end
end

PMM              = zeros(11,2);
for i=1:10
PMM(i,1)         = subM(i);       % Magnitude-sub interval

if PMM(i,1)      <= Mmin
PMM(i,2)         = 0;             % Probability Density Functions (PDFs)
else
PMM(i,2)         = fm(i);         % Probability Density Functions (PDFs)
end

end
PMM(11,1)        = Mull1;         % Mean annual rate of exceedance
PMM(11,2)        = Mull2;         % Mean annual rate of exceedance
PM               = PMM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remarks:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Supplementary scripts for print out the probability of magnitude %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pmag           = PMM(1:10,:)';
% Magname        = strcat('Output/Prob-magnitude (Fault no. ',num2str(Fnum),').txt');
% PD             = fopen(Magname,'wt');
% fprintf(PD,'%s\t%s\t' ,['Magnitude',' Probability', ]);
% fprintf(PD,'\n%f %f',[Pmag]);
% fclose(PD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
