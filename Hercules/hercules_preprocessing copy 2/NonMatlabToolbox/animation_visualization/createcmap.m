function [cmap]=createcmap(varargin)
% [cmap]=createcmap(colors)
%
% This function can be used to build your own custom colormaps. 'colors'
% can be n dim[3x1] vectors createcmap([1x3],[1x3}...) whose elements
% represent a color in RGB color space.
% Optionaly, this function can be use to build any
% colormap using main colors rgbcmyk. In image processing, w (white) can be
% used as the first color so that in the output, the background (usually
% with 0 values) appears white. In the example of rainfall map, 'wb' will
% produce a rainfall density map where the background (if its DN values are
% 0) will appear as white.
%
% Inputs:
%  colors: array ('RGB',array) of color RGB, ([1x3],[1x3]...)
%  colors: string (char) of color codes, any sequence of rgbcmywk
%  representing different colors (such as 'b' for blue) is acceptable. If a
%  gradient of white to blue is needed, colors would be 'wb'; a rainbow of
%  white+blue+red+green would be 'wbrg'.
%
% Example:
%  [cmap]=createcmap('wygbr');
% %try the output cmap:
% im=imread('cameraman.tif');
% imshow(im), colorbar
% colormap(cmap) %will use the output colormap
%
% First version: 14 Feb. 2013, modifications 14 Feb. 2015
% sohrabinia.m@gmail.com
% zu.alan.zu@gmail.com
%--------------------------------------------------------------------------

if nargin<1
    colors='wrgbcmyk';
end

if ischar(varargin{1})
    colors=varargin{1};
    ncolors=length(colors)-1;
    
    
    bins=round(255/ncolors);
    % diff1=255-bins*ncolors;
    
    vec=zeros(300,3);
    
    switch colors(1)
        case 'w'
            vec(1,:)=1;
        case 'r'
            vec(1,:)=[1 0 0];
        case 'g'
            vec(1,:)=[0 1 0];
        case 'b'
            vec(1,:)=[0 0 1];
        case 'c'
            vec(1,:)=[0 1 1];
        case 'm'
            vec(1,:)=[1 0 1];
        case 'y'
            vec(1,:)=[1 1 0];
        case 'k'
            vec(1,:)=[0 0 0];
    end
    
    
    for i=1:ncolors
        beG=(i-1)*bins+1;
        enD=i*bins+1; %beG,enD
        switch colors(i+1)
            case 'w'
                vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
                vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
                vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';%colors(i+1),beG,enD,
            case 'r'
                vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
                vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
                vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';%colors(i+1),beG,enD
            case 'g'
                vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
                vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
                vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';%colors(i+1),beG,enD
            case 'b'
                vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
                vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
                vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';%colors(i+1),beG,enD
            case 'c'
                vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
                vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
                vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';%colors(i+1),beG,enD
            case 'm'
                vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
                vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
                vec(beG:enD,3)=linspace(vec(beG,3),1,bins+1)';
            case 'y'
                vec(beG:enD,1)=linspace(vec(beG,1),1,bins+1)';
                vec(beG:enD,2)=linspace(vec(beG,2),1,bins+1)';
                vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';
            case 'k'
                vec(beG:enD,1)=linspace(vec(beG,1),0,bins+1)';
                vec(beG:enD,2)=linspace(vec(beG,2),0,bins+1)';
                vec(beG:enD,3)=linspace(vec(beG,3),0,bins+1)';
        end
    end
    cmap=vec(1:bins*ncolors,:);
    
else
    
    for kk=1:nargin
        if size(varargin{kk},2) ~= 3
         error(['Error! colors must be a variable of type char with '...
        'color-names, such as ''r'', ''g'', etc., or inputs of'...
        'at least one 1x3 array'...
        'type ''help createcmap'' for more info']);
        end
    end
    ncolors=nargin-1;
    bins=round(255/ncolors);
    vec=zeros(300,3);
    vec(1,:)=varargin{1};
    
    for i=1:ncolors
        beG=(i-1)*bins+1;
        enD=i*bins+1; %beG,enD
        Colors=varargin{i+1};
        vec(beG:enD,1)=linspace(vec(beG,1),Colors(1),bins+1)';
        vec(beG:enD,2)=linspace(vec(beG,2),Colors(2),bins+1)';
        vec(beG:enD,3)=linspace(vec(beG,3),Colors(3),bins+1)';%colors(i+1),beG,enD,
        
    end
    cmap=vec(1:bins*ncolors,:);
end
end %end of buildcmap