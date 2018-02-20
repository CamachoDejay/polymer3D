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
szWindow = 6;
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
% Finds all the position that occurs in the stack ==> some molecule are
% seen only after the Xth frame because they're out of focus before, etc...
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
zPos = repmat(z_space,size(totPos,1),1);
%% ROI extraction and Fitting
fitPar = zeros(size(totPos,1),size(imStack,3),6);
%gaussPar = zeros(size(totPos,1),size(imStack,3),3);
GraR = 4;
 im_in = double(imStack);
 
Grad = zeros(size(totPos,1),size(imStack,3));
for i=1:size(imStack,3)    
    for j=1:size(totPos,1)
        
        %Extract a roi around the localized emitter
        [roi_lims] = EmitterSim.getROI(totPos(j,1), totPos(j,2), szWindow, xSize, ySize);
        ROI = im_in(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2),i);
        ROI = imgaussfilt(ROI,2);
        %Gradient Fitting
        [x,y,e,centOut] = Localization.gradFit(ROI,GraR);
        
        %Gradient to determine focus
        [grad,~] = imgradient(ROI);
        Grad(j,i) = max(max(grad,[],2),[],1);
        
        %remove non-sense values
        if abs(x) > GraR || abs(y) > GraR
            x = NaN;
            y = NaN;
            e = 0;
        end
        %Get back the position in the original 
        xc = (round(totPos(j,1)) + x);%in px
        yc = (round(totPos(j,2)) + y);
        
        fitPar(j,i,1) = xc;
        fitPar(j,i,2) = yc;
        fitPar(j,i,3) = e;
        fitPar(j,i,4) = centOut(1).x;
        fitPar(j,i,5) = centOut(1).y;
        fitPar(j,i,6) = centOut(1).e;
    end
end
%% Forcing focal point to be z = 0 & Extracting data point
%using Fit
elipAxis = [];
zAxis    = [];
figure
hold on
for j=1:size(totPos,1)
    [out,Fit] = Misc.gauss1DFit(Grad(j,:),zPos(j,:));
    zPos(j,:) = zPos(j,:) - out(2);
    plot(zPos(j,:),fitPar(j,:,3))
    currentPar = fitPar(j,:,3);
    currentZ   = zPos(j,:);
    elipAxis   = [elipAxis, currentPar(and(and(currentZ>-1000,currentZ<1000),currentPar~=0))];
    zAxis      = [zAxis, currentZ(and(and(currentZ>-1000,currentZ<1000),currentPar~=0))];
end

hold off
%using centroid
% figure
% hold on
% for j=1:size(totPos,1)
%     plot(zPos(j,:),fitPar(j,:,6))
% end
% hold off

%%
figure
scatter(zAxis,elipAxis);



