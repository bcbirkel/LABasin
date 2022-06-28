function cmap=loadcmap(File)
%Just god and me know what this function does.
% CAMP = loadcmap(File) reads the colormap file downloaded from 
% http://soliton.vm.bytemark.co.uk/pub/cpt-city/ and creates the colormap
% for that color palete. cmap is a 255x3 matrix with the RGB colorcodes.
% Example:
% cmap=loadcmap('BlWhRe.c3g');
% Colormap=[loadcmap('tpglpom.c3g'); loadcmap('gray_lt.c3g')];
% Some downloaded colormaps:
% BlueWhiteOrangeRed, BlWhRe, cbacSpectral09, French_flag_smooth, GMT_seis,
% hx-120-120, rainbow, temperature, topography1.

FileId=fopen(File);

% If the file has not been modified here we are goin to read header:
for i=1:6
    fgetl(FileId);
end

cmap=zeros([101 3]);
n=1;
while ~feof(FileId)
   Line=fgetl(FileId);
   k=findstr('rgb',Line);
   if ~isempty(k)
       RGB(n,:)=[str2double(Line(7:9)) str2double(Line(11:13)) str2double(Line(15:17))];
       P100(n)=round(str2double(Line(20:end-2)))+1;
       n=n+1;
   end
end

RGB=RGB/255;

for k=1:length(RGB(:,1))-1
    cmap(P100(k):P100(k+1),1)=linspace(RGB(k,1),RGB(k+1,1),P100(k+1)-P100(k)+1);
    cmap(P100(k):P100(k+1),2)=linspace(RGB(k,2),RGB(k+1,2),P100(k+1)-P100(k)+1);
    cmap(P100(k):P100(k+1),3)=linspace(RGB(k,3),RGB(k+1,3),P100(k+1)-P100(k)+1);
end