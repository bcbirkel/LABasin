% distance_tau.m
% function returns the distance between two latitude and longitude points
% in degrees and kilometers
% usage:
% [deg,km] = distance([lat1 long1], [lat2 long2])
function [deg,km] = distance_tau(point1, point2)
    % check usage
    if (size(point1) ~= [1 2]) 
        disp('Size of input point1 and point2 incorrect. Should be 1x2!')
        return
    end
    if (size(point2) ~= [1 2]) 
        disp('Size of input point1 and point2 incorrect. Should be 1x2!')
        return
    end
    
    
    % get and check the inputs
    lat1 = point1(1);
    long1 = point1(2);
    lat2 = point2(1);
    long2 = point2(2);
    if (lat1 < -90 || lat1 > 90)
        disp('latitude must be between -90 and 90')
        disp('Assuming longitude and latitude reversed')
        temp=lat1;
        lat1=long1;
        long1=temp;
    end
    if (lat2 < -90 || lat2 > 90)
        disp('latitude must be between -90 and 90')
        disp('Assuming longitude and latitude reversed')
        temp=lat2;
        lat2=long2;
        long2=temp;
    end
    
    % convert to radians
    lat1rad = lat1 * pi /180;
    long1rad = long1 * pi/180;
    lat2rad = lat2 * pi/180;
    long2rad = long2 * pi/180;
    
    % compute distance in degrees
    deg = acos(sin(lat1rad) * sin(lat2rad) + cos(lat1rad) * cos(lat2rad) * cos((long2rad - long1rad))) * 180 / pi;
    
    % assuming a flat earth, compute distance in km
    km = deg * 111.19;
    
    
    return
