%preprocess_waveform_data.m
% script to read a sac file and its pole zero response then convert to ground motion
% runs in an interactive mode prompting user for file names and number of differentiations

%query user for sac file name
sacFileName=input('Please input the sac file name\n','s');

%read the sac file and store information in structure
[tmpHeader, tmpData] = load_sac(sacFileName);

%compute the nyquest frequency from the header
nyquistFreq = 1/(2*1/tmpHeader.delta);

%prepare for the fft by removing near 0 offsets
tmpData = rmean(tmpData);
tmpData = rtrend(tmpData);
tmpData = cos_taper(tmpData);

%query the user for the pole zero file name
poleZeroFileName = input('Please input the pole zero file name\n','s');

%remove the pole zero response. Assumes pole zero file returns displacement
dispData = transfer(tmpData,tmpHeader.npts, 1/tmpHeader.delta, poleZeroFileName, 0.001, nyquistFreq * 0.8 );

%query for the number of differentiations which is after assuming a displacement response
ndiff = input('Please input the number of differentiations you would like\n');

%call the diff function
for i=1:ndiff
	tmpData = diff(tmpData);
end

%let the user know the data is ready and they should copy to more stable variable
disp('Data returned in variable tmpData. Header stored in tmpHeader. Please convert these to more stable variables\n');
disp('To change the name of a variable use the syntax:\nnewVariable = oldVariable;\n');


