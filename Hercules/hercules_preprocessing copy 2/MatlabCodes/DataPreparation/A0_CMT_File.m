% Creates EarthquakeDatabase needed for later preprocessing scripts.
% CMT_Catalogue must be updated with event manually
% Modified by Brianna Birkel, May 2020

% Create mat file with source parameters from cmt file dowloaded from IRIS
clear; close all; clc;

%% Load psmeca file and parameters
%%% lon lat depth mrr mtt mpp mrt mrp mtp iexp year month day hour min
CMT_Data = load('/Users/bcbirkel/Documents/GitHub/LABasin/Hercules/hercules_preprocessing/MatlabCodes/AuxFiles/CMT_Catalogue.txt');
    
%% Assign data
for nS = 1:size(CMT_Data,1)
    % Location
    Eqk(nS,1).Lon =  CMT_Data(nS,1);
    Eqk(nS,1).Lat =  CMT_Data(nS,2);
    Eqk(nS,1).Depth =  CMT_Data(nS,3);
    Eqk(nS,1).Year =  CMT_Data(nS,11);
    Eqk(nS,1).Month =  CMT_Data(nS,12);
    Eqk(nS,1).Day =  CMT_Data(nS,13);
    Eqk(nS,1).Hour =  CMT_Data(nS,14);
    Eqk(nS,1).Min =  CMT_Data(nS,15);
    Eqk(nS,1).Sec =  CMT_Data(nS,16);
    
    
    % Moment tensor
    Eqk(nS,1).M(1,1) =  CMT_Data(nS,4);
    Eqk(nS,1).M(2,2) =  CMT_Data(nS,5);
    Eqk(nS,1).M(3,3) =  CMT_Data(nS,6);
    
    Eqk(nS,1).M(1,2) =  CMT_Data(nS,7);
    Eqk(nS,1).M(2,1) =  Eqk(nS).M(1,2);
    
    Eqk(nS,1).M(1,3) =  CMT_Data(nS,8);
    Eqk(nS,1).M(3,1) =  Eqk(nS).M(1,3);
    
    Eqk(nS,1).M(2,3) =  CMT_Data(nS,9);
    Eqk(nS,1).M(3,2) =  Eqk(nS).M(2,3);
    
    % Convert (Dyn-cm) to (Nm)
    Eqk(nS,1).M = Eqk(nS).M*(10^CMT_Data(nS,10))*(1e-7);
    
    % Estimate Julian day
    [~,S_JDay] = jl_date(Eqk(nS).Year,Eqk(nS).Month,Eqk(nS).Day);
    Eqk(nS,1).JDay = S_JDay;
end

% Save Earthquake Structure
save EarthquakeDatabase.mat Eqk;
