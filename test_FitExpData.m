% 
%The aim of this code is to receive a z-stack/ a serie of z-stack and to
%calculate the relation between the z-position and the ellipticity to serva
%as a calibration in 3D particle tracking.

% It is currently written to take .mat file (the file are prealably cropped
% as we do not consider multiplane data yet.

clear;
clc;
close all;

%% User input 
toAnalyze = 'file';
filter    = true; %use or not pre-processing Gaussian filter
FWHM_nm   = 350;%in nm
pxSize    = 105;%in nm
szWindow  = 6;
zSpacing  = 50; %in nm

setupInfo.FWHM     = FWHM_nm;
setupInfo.pxSize   = pxSize;
setupInfo.szWindow = szWindow;
setupInfo.zCalibration = [-674.06 677.81 -77.022];
%% Loading of the data
switch toAnalyze
    case 'file'
        [fileName,folder]=uigetfile({'*.mat'},'Select a file to Open');
        [~,fileTif,~]=fileparts(fileName);
        [~,fileOME,~]=fileparts(fileTif);

        path2data = [folder fileName];

        tmp=load(path2data);
        name=fields(tmp);
        imStack=tmp.(name{1});

        %% Z-Calibration
        %[zAxis,ellipAxis] = zCalibration(setupInfo,imStack,filter);
        [Loc] = Fitting_ExpData(setupInfo,imStack);

%         figure
%         scatter(zAxis,ellipAxis);
    case 'folder'
        mainFolderName = uigetdir;
        assert(ischar(mainFolderName),'User canceled the selection of file, excecution aborted');

        %Extract the part of the folder that is a tif file
        Folder_Content = dir(mainFolderName);
        index2Images      = contains({Folder_Content.name},'.mat');
        images2Analyze = Folder_Content(index2Images);
        
        Results(size(images2Analyze,1)).fileName = [];
        Results(size(images2Analyze,1)).xc = [];
        Results(size(images2Analyze,1)).yc = [];
        Results(size(images2Analyze,1)).zc = [];
        Results(size(images2Analyze,1)).label = [];
        h = waitbar(0,'Localization...');
        for i=1:size(images2Analyze,1)
            
            file2Load = sprintf('%s%s%s',images2Analyze(i).folder,'\',images2Analyze(i).name);
            tmp = load(file2Load);
            name=cellstr(fields(tmp));
            imStack=tmp.(name{1});
            [Loc] = Fitting_ExpData(setupInfo,imStack);
            Results(i).fileName = images2Analyze(i).name;
            Results(i).xc = Loc.xc;
            Results(i).yc = Loc.yc;
            Results(i).zc = Loc.zc;
            Results(i).label = Loc.label;
            waitbar(i/size(images2Analyze,1),h);
        end
        close(h);
end


