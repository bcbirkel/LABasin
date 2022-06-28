function CopyRequiredFiles
% This function allows the user to select a matlab function and
% destination folder. It will then copy any required functions to the
% destination folder.

clear,clc

Cancel = 0;
% SELECT MASTER MATLAB FILE
[Name, Location] = uigetfile({'*.m'},'Select master matlab file');
Source = fullfile(Location,Name);


if Location ~= 0
    % SELECT DESTINATION FOLDER
    Destination = uigetdir(Location,'Select folder to copy to');

    if Destination ~= 0
        % GET REQUIRED FILES AND LOCATIONS
        disp('Acquiring required files')
        CurrentFolder = pwd;
        [pathstr,~,~] = fileparts(Source);
        cd(pathstr);
        Functions = matlab.codetools.requiredFilesAndProducts(Source);
        cd(CurrentFolder);
        disp('Files acquired')
        NF = length(Functions);

        NC = zeros(1,NF);
        % CHECK IF FILES EXIST
        for i = 1:NF
            [~,name,ext] = fileparts(char(Functions(i)));
            if exist([Destination,'\',name,ext],'file') ~= 0
                button = questdlg(['Overwrite file: ',Destination,'\',name,ext],'FWI','Yes','Yes to all remaining','No','No');

                if strcmp('Yes',button) == 1
                    delete([Destination,'\',name,ext])
                    NC(i) = 2;

                elseif strcmp('Yes to all remaining',button) == 1
                    for ii = i:NF
                        [~,name2,ext2] = fileparts(char(Functions(ii)));
                        if exist([Destination,'\',name2,ext2],'file') ~= 0
                            delete([Destination,'\',name2,ext2])
                            NC(ii) = 2;
                        else
                            NC(ii) = 1;
                        end
                    end
                    break

                elseif strcmp('No',button) == 1
                    NC(i) = 0;
                    
                else
                    Cancel = 1;
                    break
                end
            else
                NC(i) = 1;
            end
        end
        % COPY FILES
        if Cancel == 0
            for i = 1:NF
                [~,name,ext] = fileparts(char(Functions(i)));
                if NC(i) == 0
                    disp(['Not Copied: ',Destination,'\',name,ext])
                end
                if NC(i) == 1
                    copyfile(char(Functions(i)),Destination)
                    disp(['Copied: ',Destination,'\',name,ext])
                end
                if NC(i) == 2
                    copyfile(char(Functions(i)),Destination)
                    disp(['Overwritten: ',Destination,'\',name,ext])
                end
            end
            disp('Copy complete')
        else
            disp('Copy canceled')
        end
    else 
        disp('Copy canceled')
    end
else 
    disp('Copy canceled')
end
