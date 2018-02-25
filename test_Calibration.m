%%
%The aim of this code is to receive a z-stack/ a serie of z-stack and to
%calculate the relation between the z-position and the ellipticity to serva
%as a calibration in 3D particle tracking.

% It is currently written to take .mat file (the file are prealably cropped
% as we do not consider multiplane data yet.

clear;
clc;
close all;

%% User input
toAnalyze = 'folder';
filter    = false; %use or not pre-processing Gaussian filter
FWHM_nm   = 350;%in nm
pxSize    = 105;%in nm
szWindow  = 6;
zSpacing  = 50; %in nm

setupInfo.FWHM     = FWHM_nm;
setupInfo.pxSize   = pxSize;
setupInfo.szWindow = szWindow;
setupInfo.zSpacing = zSpacing;
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
        [zAxis,ellipAxis] = zCalibration(setupInfo,imStack,filter);
        
        figure
        scatter(zAxis,ellipAxis);
    case 'folder'
        mainFolderName = uigetdir;
        assert(ischar(mainFolderName),'User canceled the selection of file, excecution aborted');
        
        %Extract the part of the folder that is a tif file
        Folder_Content = dir(mainFolderName);
        index2Images      = contains({Folder_Content.name},'.mat');
        images2Analyze = Folder_Content(index2Images);
        
        zAxis(size(images2Analyze,1)).cam1 = [];
        zAxis(size(images2Analyze,1)).cam2 = [];
        zAxis(size(images2Analyze,1)).cam1All = [];
        zAxis(size(images2Analyze,1)).cam2All = [];
        zAxis(size(images2Analyze,1)).all = [];
        ellipAxis(size(images2Analyze,1)).cam1 = [];
        ellipAxis(size(images2Analyze,1)).cam2 = [];
        ellipAxis(size(images2Analyze,1)).cam1All = [];
        ellipAxis(size(images2Analyze,1)).cam2All = [];
        ellipAxis(size(images2Analyze,1)).all = [];
        h = waitbar(0,'Z-Calibration...');
        for i=1:size(images2Analyze,1)
            
            file2Load = sprintf('%s%s%s',images2Analyze(i).folder,'\',images2Analyze(i).name);
            tmp = load(file2Load);
            name=cellstr(fields(tmp));
            imStack=tmp.(name{1});
            [z,ellip] = zCalibration(setupInfo,imStack,filter);
            
            if contains(file2Load,'Cam1')
                zAxis(i).cam1     = z;
                ellipAxis(i).cam1 = ellip;
                zAxis(1).cam1All = [[zAxis(1).cam1All] z];
                ellipAxis(1).cam1All = [[ellipAxis(1).cam1All] ellip];
            elseif contains(file2Load,'Cam2')
                zAxis(i).cam2     = z;
                ellipAxis(i).cam2 = ellip;
                zAxis(1).cam2All = [[zAxis(1).cam2All] z];
                ellipAxis(1).cam2All = [[ellipAxis(1).cam2All] ellip];
            end
            zAxis(1).all = [[zAxis(1).all] z];
            ellipAxis(1).all = [[ellipAxis(1).all] ellip];
            waitbar(i/size(images2Analyze,1),h);
        end
        %Sorting
        [zAxis(1).all,B] = sort([zAxis(1).all]);
        ellipAxis(1).all = [ellipAxis(1).all(B)];
        
        [zAxis(1).cam1All,C] = sort([zAxis(1).cam1All]);
        ellipAxis(1).cam1All = [ellipAxis(1).cam1All(C)];
        
        [zAxis(1).cam2All,D] = sort([zAxis(1).cam2All]);
        ellipAxis(1).cam2All = [ellipAxis(1).cam2All(D)];
        
        figure
        scatter(zAxis(1).all,ellipAxis(1).all);
        
end
close(h);



