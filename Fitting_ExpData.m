function [locPos] = Fitting_ExpData(setupInfo,imStack)

FWHM_nm = setupInfo.FWHM;
pxSize  = setupInfo.pxSize;
szWindow = setupInfo.szWindow;
zCalibration1 = setupInfo.zCalibration1;
zCalibration2 = setupInfo.zCalibration2;
xSize = size(imStack,2);
ySize = size(imStack,1);

%% Localization
% Finds all the position that occurs in the stack ==> some molecule are
% seen only after the Xth frame because they're out of focus before, etc...

for i = 1:size(imStack,3)
    im_in = double(imStack(:,:,i));
    delta = 4;
    FWHM_pix = FWHM_nm / pxSize; %[pix]
    % for GLRT
    chi2 = 24;
    [ pos, ~] = Localization.smDetection(im_in, delta, FWHM_pix, chi2 );
    if i==1
        totPos = pos;
        cmpPos = pos;
        label  = 1:size(pos,1);
        label  = label(:);
    else
        for j=1:size(pos,1)
            testX = and(cmpPos(:,1) < pos(j,1)+8,cmpPos(:,1) > pos(j,1)-8);
            testY = and(cmpPos(:,2) < pos(j,2)+8,cmpPos(:,2) > pos(j,2)-8);
            idxX  = find(testX == 1);
            idxY  = find(testY == 1);
            add2Pos = true;
            if(~isempty(idxX)) && ~isempty(idxY)
                for k = 1:size(idxX,1)
                    testidxY = idxY == idxX(k);
                    if(size(testidxY(testidxY==1),1)>1)
                        add2Pos = false;
                    elseif (size(testidxY(testidxY==1),1)==1)
                        testidx = testidxY;
                    end
                end
            else
                add2Pos = false;
            end
            if(add2Pos)
                totPos(end+1,:) = pos(j,:);
                label(end+1,:)  = label(idxY(testidx));
 
            end
            
        end
    end
end

%% ROI extraction and Fitting
fitPar = zeros(size(imStack,3),max(label),8);
%gaussPar = zeros(size(totPos,1),size(imStack,3),3);
GraR = 4;
im_in = double(imStack);
    Frame = 1;
    for j=1:size(label,1)%max num of molecule found
        if j~=size(label,1) && j>1
            if and(label(j-1)>ceil(max(label)/2),label(j)<=floor(max(label)/2))
                Frame = Frame+1;
            end
%         else
%             Frame = size(label,1);
        end
           %dbstop in Fitting_ExpData at 70 if Frame==11
        %Extract a roi around the localized emitter
        [roi_lims] = EmitterSim.getROI(totPos(j,1), totPos(j,2), szWindow, xSize, ySize);
        ROI = im_in(roi_lims(3):roi_lims(4),roi_lims(1):roi_lims(2),Frame);
        
        %Gradient Fitting
        if size(ROI,1)~= szWindow*2+1 || size(ROI,2)~= szWindow*2+1
            x = NaN;
            y = NaN;
            e = NaN;
            centOut.x = NaN;
            centOut.y = NaN;
            centOut.e = NaN;
        else
            [x1,y1,~,~] = Localization.gradFit(ROI,GraR);
            ROI = imgaussfilt(ROI,2);
            [~,~,e1,~] = Localization.gradFit(ROI,GraR);
            [x,y,e] = Localization.phasor(ROI);
        end

        %remove non-sense values
        if abs(x) > GraR || abs(y) > GraR || e <=0
            x = NaN;
            y = NaN;
            e = NaN;
        end
        %Get back the position in the original
        x1c =(round(totPos(j,1)) + x1);
        y1c =(round(totPos(j,2)) + y1);
        xc = (round(totPos(j,1)) + x);%in px
        yc = (round(totPos(j,2)) + y);
        
        
        fitPar(Frame,label(j),1) = xc;
        fitPar(Frame,label(j),2) = yc;
        fitPar(Frame,label(j),5) = x1c;
        fitPar(Frame,label(j),6) = y1c;
        fitPar(Frame,label(j),4) = e;
        fitPar(Frame,label(j),8) = e1;
      
        
        %Calculating z:
        fitPar(Frame,label(j),3) = zCalibration2(3)*e^2+zCalibration2(2)*e+zCalibration2(1);
        fitPar(Frame,label(j),7) = zCalibration1(3)*e1^2+zCalibration1(2)*e1+zCalibration1(1);
        
    end

locPos.xc = fitPar(:,:,1);
locPos.yc = fitPar(:,:,2);
locPos.zc = fitPar(:,:,3);
locPos.e  = fitPar(:,:,4);
locPos.x1c = fitPar(:,:,5);
locPos.y1c = fitPar(:,:,6);
locPos.z1c = fitPar(:,:,7);
locPos.e1  = fitPar(:,:,8);
locPos.label = label;


end