% this script is a work in progress, the plan is to fix here the quick
% callibration of the channels. Basically I want to be able to lead the
% movies and built a dataset of [xDim,yDim,channel,frame] where all
% channels represent the same region of the image to a pixel resolution.
clear
close all
clc

fPath = '/Users/rafa/Documents/MATLAB/data/multi-plane/Beads - MC_4';
fName = 'Beads - MC_4_MMStack_Pos0.ome.tif';

fPath = [fPath filesep fName];

% Calculate calibration
[cal] = mpSetup.cali.calculate(fPath, false);

% load and calibrate, when applied to the calibration data then we should
% be able to demonstrate that it works
[data] = mpSetup.loadAndCal( fPath, cal );

% figure to show that it works
f = round(mean(cat(1,cal.inFocus.frame)));
xl = [97,206];
yl = [116,229];
for i =1:8
    subplot(2,4,i)
    imagesc(data(:,:,i,f))
    axis image
%     xlim(xl)
%     ylim(yl)
%     hold on
    grid on
    
    title(['Ordered channel ' num2str(i)])
    a = gca;
    a.XTickLabel = [];
    a.YTickLabel = [];
end



