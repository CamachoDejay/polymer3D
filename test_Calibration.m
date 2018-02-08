%% 
%The aim of this code is to receive a z-stack/ a serie of z-stack and to
%calculate the relation between the z-position and the ellipticity to serva
%as a calibration in 3D particle tracking.

% It is currently written to take .mat file (the file are prealably cropped
% as we do no consider multiplane data yet.

clear;
clc;
close all;

%% User input 
FWHM_nm = 350;%in nm
pxSize  = 105;%in nm
z_spacing = 50; %in nm
szWindow =6;
%% Loading of the data
[fileName,folder]=uigetfile({'*.mat'},'Select a file to Open');
[~,fileTif,~]=fileparts(fileName);
[~,fileOME,~]=fileparts(fileTif);

path2data = [folder fileName];

tmp=load(path2data);
name=fields(tmp);
imStack=tmp.(name{1});

xSize = size(imStack,2);
ySize = size(imStack,1);
z_space = linspace(0,size(imStack,3)*z_spacing,size(imStack,3));
%% Localization
% Finds all the position that occurs in the stack
totPos = [0,0];
for i = 1:size(imStack,3)
    im_in = double(imStack(:,:,i));
     delta = 4;
     FWHM_pix = FWHM_nm / pxSize; %[pix]
     % for GLRT
     chi2 = 24;
[ pos, inten ] = Localization.smDetection(im_in, delta, FWHM_pix, chi2 );

    for j=1:size(pos,1)
        testX = and(totPos(:,1) < pos(j,1)+4,totPos(:,1) > pos(j,1)-4);
        testY = and(totPos(:,2) < pos(j,2)+4,totPos(:,2) > pos(j,2)-4);
        idxX  = find(testX == 1);
        idxY  = find(testY == 1);
        add2Pos = true;
        if(~isempty(idxX))
            for k = 1:size(idxX)
            testidxY = idxY == idxX(k);
            if(size(testidxY==1)>=1)
                add2Pos = false;
            end
            end
       end
       if(add2Pos)
           totPos(end+1,:) = pos(j,:);
       end
    end
end
totPos(1,:)=[];

%% ROI extraction and Fitting
fitPar = zeros(size(totPos,1),size(imStack,3),6);
%gaussPar = zeros(size(totPos,1),size(imStack,3),3);
GraR = 4;
 im_in = double(imStack);
for i=1:size(imStack,3)    
    for j=1:size(totPos,1)
        
        %Extract a roi around the localized emitter
        [roi_lims] = EmitterSim.getROI(totPos(j,1), totPos(j,2), szWindow, xSize, ySize);
        ROI = im_in(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2),i);
        
        [x,y,e,centOut] = Localization.gradFit(ROI,GraR);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %num and den used for calculating e have non-sense values (10^21) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Extracting ellipticity from 1D Gaussian Fit
%         XData = max(ROI,[],1);
%         YData = max(ROI,[],2)';
%         rawDataX.N = XData;
%         rawDataX.binCenter = 1:1:length(XData);
%         paramsX = [FWHM_pix/2.3,round(length(XData)/2),max(XData),min(XData)]; 
%         paramsY = [FWHM_pix/2.3,round(length(YData)/2),max(YData),min(YData)]; 
%         rawDataY.N = YData;
%         rawDataY.binCenter = 1:1:length(YData);
        
%         funX = @(x) GradientFit.fun2minGauss(rawDataX,x,false);
%         funY = @(x) GradientFit.fun2minGauss(rawDataY,x,false);
%         % then we can do:
%         [outX, RMSDX] = fminsearch(funX,paramsX);
%         [outY,RMSDY]  = fminsearch(funY,paramsY);
%         
%         gaussPar(j,i,1) = outX(1);
%         gaussPar(j,i,2) = outY(1);
%         gaussPar(j,i,3) = outX(1)/outY(1);
        
        if abs(x) > GraR || abs(y) > GraR
            x = NaN;
            y = NaN;
            e = 0;
        end
        
        xc = (round(totPos(j,1)) + x);%in px
        yc = (round(totPos(j,2)) + y);
        
        fitPar(j,i,1) = xc;
        fitPar(j,i,2) = yc;
        fitPar(j,i,3) = e;
        fitPar(j,i,4) = centOut.x;
        fitPar(j,i,5) = centOut.y;
        fitPar(j,i,6) = centOut.e;
    end
end
