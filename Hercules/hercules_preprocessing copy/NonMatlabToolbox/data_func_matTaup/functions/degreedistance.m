function [dist] = degreedistance(ll)
%degreedistance.m
%Function to find the angular distance between two [lat,long] sites in
%degrees and distance in kilometers
%Usage: degreedistance([lat1,long1;lat2,long2]}
xy1=lltoxy([ll(1,1),ll(1,2)]);
xy2=lltoxy([ll(2,1),ll(2,2)]);

distance=acos(dot(xy1,xy2))*180/pi;
kmdistance=distance * 111.19;

dist=[distance;kmdistance];