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
FWHM_nm = 200;
FWHM_pix = FWHM_nm / pxSize;
% for GLRT
chi2 = 24;

tPos{8,1}=[];

figure(1)
clf

for i = 1:8
    im = data(:,:,i);
    im = double(im);
    [ pos, inten ] = Localization.smDetection( im, delta, FWHM_pix, chi2 );
    tPos{i} = [pos,inten];
    subplot(2,4,i)
    imagesc(im)
    axis image
    hold on
    plot(pos(:,1),pos(:,2),'xr')
    hold off
end

% make a common list of sm detections? should I have a test for seeing a
% molecule in at least 3 (or X) planes? once I have the common list I have
% to cropt the ROIs [dx, dy, 8] and do the fine localization. We will for
% the moment pick the best as I do not have a super-resolved registration
% matrix between all channels. This should be generated in order to use
% info from multiple planes to increase fit accuracy. it is from z tacks of
% beads that we can create such a registration.



