function [output] = zCalibration(setupInfo,imStack,filter,method)

FWHM_nm = setupInfo.FWHM;
pxSize  = setupInfo.pxSize;
szWindow = setupInfo.szWindow;
zSpacing = setupInfo.zSpacing; %in nm

xSize = size(imStack,2);
ySize = size(imStack,1);
z_space = linspace(0,size(imStack,3)*zSpacing,size(imStack,3));
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
    [ pos, ~] = Localization.smDetection(im_in, delta, FWHM_pix, chi2 );
    
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
        
        
        
        %Gradient Fitting
        if size(ROI,1)~= szWindow*2+1 || size(ROI,2)~= szWindow*2+1
            x = NaN;
            y = NaN;
            e = NaN;
            centOut.x = NaN;
            centOut.y = NaN;
            centOut.e = NaN;
        else
            switch method
                case 'gradient'
                    [x,y,e,centOut] = Localization.gradFit(ROI,GraR);
                    if filter
                        ROI = imgaussfilt(ROI,2);
                        [~,~,e,~] = Localization.gradFit(ROI,GraR);
                    end
                otherwise
                    [x,y,e] = Localization.phasor(ROI);
                    centOut.x = NaN;
                    centOut.y = NaN;
                    centOut.e = NaN;
            end
            
        end
        %Gradient to determine focus
        [grad,~] = imgradient(ROI);
        Grad(j,i) = max(max(grad,[],2),[],1);
        
        %remove non-sense values
        if abs(x) > GraR || abs(y) > GraR || e <=0
            x = NaN;
            y = NaN;
            e = NaN;
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

figure
hold on
for j=1:size(totPos,1)
    [out,fit] = Misc.gauss1DFit(Grad(j,:),zPos(j,:));
    plot(Grad(j,:),'--')
    plot(fit)
    zPos(j,:) = zPos(j,:) - out(2);
    currentPar = fitPar(j,:,3);
    currentZ   = zPos(j,:);
    currentX   = fitPar(j,:,1);
    currentY   = fitPar(j,:,2);
    ellip{j}   = currentPar(and(and(currentZ>-1000,currentZ<1000),currentPar~=0));
    zcell{j}   = currentZ(and(and(currentZ>-1000,currentZ<1000),currentPar~=0));
    xcell{j}   = currentX(and(and(currentZ>-1000,currentZ<1000),currentPar~=0));
    ycell{j}   = currentY(and(and(currentZ>-1000,currentZ<1000),currentPar~=0));
end

output.x = xcell;
output.y = ycell;
output.z = zcell;
output.ellip = ellip;
end