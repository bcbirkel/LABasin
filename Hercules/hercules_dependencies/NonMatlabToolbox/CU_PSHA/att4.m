function PGASD = att4(M,R,Ft)
%-------------------------------------------------------------------------%
%                                                                         %
%        Strong Ground-Motion Attenuation Relationship (PGA)              %
%                                                                         %
% Reference: Boore et al. (1993) & Boore et al. (1997)                    %
%-------------------------------------------------------------------------%

%(Step 1) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Set Coefficient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for randomly-oriented horizontal component (or geometrical mean) 
b1         = -0.105;
b2         = 0.229; 
b3         = 0; 
b4         = 0; 
b5         = -0.778; 
b6         = 0.162; 
b7         = 0.251; 
h          = 5.57;
SD         = 0.230;
% Use three site categories:
% Class A Vs,30 > 750 ms1, some categorised using measured shear-wave velocity, most estimated )
    GB     = 0;
    GC     = 0;
% Class B 360 < Vs,30  750 ms1, some categorised using measured shear-wave velocity, most estimated ) 
%   GB     = 1,GC = 0, 118 records.
% Class C 180 < Vs,30  360 ms1,some categorised using measured shear-wave velocity, most estimated ) 
%   GB     = 0,GC = 1, 105 records.
% where Vs,30 is average shear-wave velocity to 30m.

%(Step 2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculate mean PGA and Standard Deviation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R1=sqrt((R^2)+(h^2));
PGA            = log((10^(b1+(b2*(M-6))+(b3*((M-6)^2))+(b4*R1)+(b5*log10(R1))+(b6*GB)+(b7*GC)))*981); % in unit "ln(gal)" (1 gal=0.01 m/s2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PGASD          = [PGA SD];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%