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
toAnalyze = 'folder';
FWHM_nm   = 350;%in nm
pxSize    = 95;%in nm
szWindow  = 6;
zSpacing  = 50; %in nm

setupInfo.FWHM     = FWHM_nm;
setupInfo.pxSize   = pxSize;
setupInfo.szWindow = szWindow;
setupInfo.zCalibration1 = [+1.077 0.00131];
setupInfo.zCalibration2 = [-674.06 677.81 -77.022];
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

        %% Fitting
        %[zAxis,ellipAxis] = zCalibration(setupInfo,imStack,filter);
        [Loc] = Fitting_ExpData(setupInfo,imStack);
        Loc.xc = (Loc.xc-Loc.xc(1))*pxSize;
        Loc.yc = (Loc.yc-Loc.yc(1))*pxSize;
        Loc.zc = (Loc.zc-Loc.zc(1));
        Loc.x1c = (Loc.x1c-Loc.x1c(1))*pxSize;
        Loc.y1c = (Loc.y1c-Loc.y1c(1))*pxSize;
        Loc.z1c = (Loc.z1c-Loc.z1c(1));
        
        figure
        plot(Loc.xc)
        hold on
        plot(Loc.yc)
        plot(Loc.x1c)
        plot(Loc.y1c)
        title('X-Y movement of beads with stage')
        xlabel('Frames')
        plot(Loc.zc)
        plot(Loc.z1c)
        legend('XPhasor', 'YPhasor', 'XGrad','YGrad','ZPhasor','ZGrad');
        ylabel('Distance (nm)')
        
        
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
        Results(size(images2Analyze,1)).e = [];
        Results(size(images2Analyze,1)).xcGrad = [];
        Results(size(images2Analyze,1)).ycGrad = [];
        Results(size(images2Analyze,1)).zcGrad = [];
        Results(size(images2Analyze,1)).e1 = [];
        Results(size(images2Analyze,1)).label = [];
        h = waitbar(0,'Localization...');
        for i=1:size(images2Analyze,1)
            
            file2Load = sprintf('%s%s%s',images2Analyze(i).folder,'\',images2Analyze(i).name);
            tmp = load(file2Load);
            name=cellstr(fields(tmp));
            imStack=tmp.(name{1});
            [Loc] = Fitting_ExpData(setupInfo,imStack);
            Results(i).fileName = images2Analyze(i).name;
            for j = 1 : size(Loc.xc,2)
            
            Results(i).xc(:,j) = (Loc.xc(:,j)-Loc.xc(1,j))*pxSize;
            Results(i).yc(:,j) = (Loc.yc(:,j)-Loc.yc(1,j))*pxSize;
            Results(i).zc(:,j) = (Loc.zc(:,j)-Loc.zc(1,j));
            Results(i).e(:,j)  = Loc.e(:,j);
            Results(i).xcGrad(:,j) = (Loc.x1c(:,j)-Loc.x1c(1,j))*pxSize;
            Results(i).ycGrad(:,j) = (Loc.y1c(:,j)-Loc.y1c(1,j))*pxSize;
            Results(i).zcGrad(:,j) = (Loc.z1c(:,j)-Loc.z1c(1,j));
            Results(i).e1(:,j)  = Loc.e1(:,j);
            end
            Results(i).label = Loc.label;
            waitbar(i/size(images2Analyze,1),h);
        end
        close(h);
end


