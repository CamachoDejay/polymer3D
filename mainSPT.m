% this script is a work in progress, the plan is to fix here the quick
% callibration of the channels. Basically I want to be able to lead the
% movies and built a dataset of [xDim,yDim,channel,frame] where all
% channels represent the same region of the image to a pixel resolution.
clear
close all
clc

%% path to the callibration
fPath = '/Users/rafa/Documents/MATLAB/data/Boris/180322-Boris-Calmultiplane/BeadsCalibrationZStack_2';
fName = 'BeadsCalibrationZStack_2_MMStack_Pos0.ome.tif';

fPath = [fPath filesep fName];

% Calculate calibration
[cal] = mpSetup.cali.calculate(fPath, false);

%%
% load and calibrate, when applied to the calibration data then we should
% be able to demonstrate that it works

fPath = '/Users/rafa/Documents/MATLAB/data/Boris/180322-Boris-Calmultiplane/TL-OD2-200msExposureR2_1';
fName = 'TL-OD2-200msExposureR2_1_MMStack_Pos0.ome.tif';
fPath = [fPath filesep fName];
[data, frameInfo, movInfo] = mpSetup.loadAndCal( fPath, cal, 1);

% figure to show that it works
% f = 1
% for i =1:8
%     subplot(2,4,i)
%     imagesc(data(:,:,i,f))
%     axis image
%     grid on
%     
%     title(['Ordered channel ' num2str(i)])
%     a = gca;
%     a.XTickLabel = [];
%     a.YTickLabel = [];
% end

%%
delta = 4;
pxSize = 100;
FWHM_nm = 300;
FWHM_pix = FWHM_nm / pxSize;
% for GLRT
chi2 = 50;

% generate a list of all localized molecules per im plane, together with a
% common list for all im planes
tLoc{8,1}=[];
rTh = 4;
consLoc = [];

for i = 1:8
    
    im = data(:,:,i);
    im = double(im);
    [ pos, inten ] = Localization.smDetection( im, delta, FWHM_pix, chi2 );
    tLoc{i} = [pos,inten];
    
    if isempty(consLoc)
        consLoc = pos;
    else
        if ~isempty (pos)
            [ consLoc ] = Localization.consolidatePos( consLoc, pos, rTh );
        end
        
    end
    
end
%%
% build ROIs
ROIcent = round(consLoc);
ROIcent(end,:) = [457,670];
molRad = 5;
nLoc   = size(ROIcent,1);
imSize = [size(data,1),size(data,2)];

% test for molecules to close to the im edges
tVal = (ROIcent(:,1) - molRad);
test = tVal < 1;
ROIcent(test,1) = molRad +1;

tVal = (ROIcent(:,2) - molRad);
test = tVal < 1;
ROIcent(test,2) = molRad +1;

tVal = (ROIcent(:,1) + molRad);
test = tVal > imSize(2);
ROIcent(test,1) = imSize(2) - molRad;

tVal = (ROIcent(:,2) + molRad);
test = tVal > imSize(1);
ROIcent(test,2) = imSize(1) - molRad;


ROIs = zeros(nLoc,4);
ROIs(:,1) = round(ROIcent(:,1))-molRad;
ROIs(:,2) = round(ROIcent(:,1))+molRad;
ROIs(:,3) = round(ROIcent(:,2))-molRad;
ROIs(:,4) = round(ROIcent(:,2))+molRad;

imLims = repmat([1 imSize(2) 1 imSize(1)],nLoc,1);
test = [ROIs(:,1)>=imLims(:,1), ROIs(:,2)<=imLims(:,2),...
        ROIs(:,3)>=imLims(:,3), ROIs(:,4)<=imLims(:,4)];
    
assert(all(test(:)), 'Problems, unexpected issues during ROI definition')


%%

% figure to check
figure(1)
clf
for i = 1:8
    im = data(:,:,i);
    pos = tLoc{i};

    subplot(2,4,i)
    imagesc(im)
    axis image
    hold on
    if ~isempty(pos)
        plot(pos(:,1),pos(:,2),'og','markersize',20)
    end
    plot(consLoc(:,1),consLoc(:,2),'xr','markersize',5) 
%     for j=1 : nLoc
%         plot(ROI(1:2,j),ROI(3:4,j),'-b')
%     end
    hold off
%     colormap hot
%     title(['Im Channel: ' num2str(i)])
%     a = gca;
%     a.FontSize = 14;
    
end

% make a common list of sm detections? should I have a test for seeing a
% molecule in at least 3 (or X) planes? once I have the common list I have
% to cropt the ROIs [dx, dy, 8] and do the fine localization. We will for
% the moment pick the best as I do not have a super-resolved registration
% matrix between all channels. This should be generated in order to use
% info from multiple planes to increase fit accuracy. it is from z tacks of
% beads that we can create such a registration.



