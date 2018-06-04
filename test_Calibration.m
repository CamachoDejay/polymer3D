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
method    = 'phasor';
filter    = false; %use or not pre-processing Gaussian filter
FWHM_nm   = 350;%in nm
pxSize    = 95;%in nm
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
        
        tmp     = load(path2data);
        name    = fields(tmp);
        imStack = tmp.(name{1});
        
        %% Z-Calibration
        [output]  = zCalibration(setupInfo,imStack,filter,method);
        zAxis =[];
        ellipAxis = [];
        for i = 1: size(output.z,2)
            zAxis = [zAxis [output.z{i}]];
            ellipAxis = [ellipAxis [output.ellip{i}]];
        end
        [zAxis,Ind] = sort(zAxis);
        ellipAxis = ellipAxis(Ind);
        
        figure
        scatter(zAxis,ellipAxis);
    case 'folder'
        mainFolderName = uigetdir;
        assert(ischar(mainFolderName),'User canceled the selection of file, excecution aborted');
        mkdir(mainFolderName,'Results');
        newdir         = sprintf('%s%s',mainFolderName,'\Results\');

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

        xAxis(size(images2Analyze,1)).cam1 = [];
        xAxis(size(images2Analyze,1)).cam2 = [];
        xAxis(size(images2Analyze,1)).cam1All = [];
        xAxis(size(images2Analyze,1)).cam2All = [];
        xAxis(size(images2Analyze,1)).all = [];

        yAxis(size(images2Analyze,1)).cam1 = [];
        yAxis(size(images2Analyze,1)).cam2 = [];
        yAxis(size(images2Analyze,1)).cam1All = [];
        yAxis(size(images2Analyze,1)).cam2All = [];
        yAxis(size(images2Analyze,1)).all = [];
        
        h = waitbar(0,'Z-Calibration...');
        for i=1:size(images2Analyze,1)
            
            file2Load = sprintf('%s%s%s',images2Analyze(i).folder,'\',images2Analyze(i).name);
            tmp = load(file2Load);
            name=cellstr(fields(tmp));
            imStack=tmp.(name{1});
            
            [output] = zCalibration(setupInfo,imStack,filter,method);
            
            if contains(file2Load,'Cam1')
                zAxis(i).cam1     = output.z;
                xAxis(i).cam1     = output.x;
                yAxis(i).cam1     = output.y;
                ellipAxis(i).cam1 = output.ellip;
                
                zAxis(1).cam1All = [[zAxis(1).cam1All] output.z];
                xAxis(1).cam1All = [[xAxis(1).cam1All] output.x];
                yAxis(1).cam1All = [[yAxis(1).cam1All] output.y];
                ellipAxis(1).cam1All = [[ellipAxis(1).cam1All] output.ellip];
                
            elseif contains(file2Load,'Cam2')
                zAxis(i).cam2     = output.z;
                xAxis(i).cam2     = output.x;
                yAxis(i).cam2     = output.y;
                ellipAxis(i).cam2 = output.ellip;
                
                zAxis(1).cam2All = [[zAxis(1).cam2All] output.z];
                xAxis(1).cam2All = [[xAxis(1).cam2All] output.x];
                yAxis(1).cam2All = [[yAxis(1).cam2All] output.y];
                ellipAxis(1).cam2All = [[ellipAxis(1).cam2All] output.ellip];
            end
            zAxis(1).all = [[zAxis(1).all] output.z];
            xAxis(1).all = [[xAxis(1).all] output.x];
            yAxis(1).all = [[yAxis(1).all] output.y];
            ellipAxis(1).all = [[ellipAxis(1).all] output.ellip];
            waitbar(i/size(images2Analyze,1),h);
        end
        close(h);
        
        %Sorting
%         [zAxis(1).all,B] = sort([zAxis(1).all]);
%         ellipAxis(1).all = [ellipAxis(1).all(B)];
%         
%         [zAxis(1).cam1All,C] = sort([zAxis(1).cam1All]);
%         ellipAxis(1).cam1All = [ellipAxis(1).cam1All(C)];
%         
%         [zAxis(1).cam2All,D] = sort([zAxis(1).cam2All]);
%         ellipAxis(1).cam2All = [ellipAxis(1).cam2All(D)];
%         
%         figure
%         scatter(zAxis(1).all,ellipAxis(1).all);
        
        fileXAxis = sprintf('%sxAxis.mat',newdir);
        fileYAxis = sprintf('%syAxis.mat',newdir);
        fileZAxis = sprintf('%szAxis.mat',newdir);
        fileEllipAxis = sprintf('%sellipAxis.mat',newdir);
        
        save(fileXAxis,'xAxis');
        save(fileYAxis,'yAxis');
        save(fileZAxis,'zAxis');
        save(fileEllipAxis,'ellipAxis');
        
end




