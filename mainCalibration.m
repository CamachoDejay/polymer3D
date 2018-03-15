% this script is a work in progress, the plan is to fix here the quick
% callibration of the channels. Basically I want to be able to lead the
% movies and built a dataset of [xDim,yDim,channel,frame] where all
% channels represent the same region of the image to a pixel resolution.
clear
close all
clc

fPath = 'D:\Documents\Unif\PhD\2018-Data\03 - March\13\ZstackCalibration-100nmBeads-640nm_1';
fName = 'ZstackCalibration-100nmBeads-640nm_1_MMStack_Pos0.ome.tif';

fPath = [fPath filesep fName];

% Calculate calibration
[cal] = mpSetup.cali.calculate(fPath, false);

%TODO: Improve the calculation of distance between the plane, Move it to
%calculate ?
zFocus = zeros(1,size(cal.focusMet,2));
zFocus2 = zeros(1,size(cal.focusMet,2));
for k=1:size(cal.focusMet,2)
    [out,Fit] = Misc.gauss1DFit(cal.focusMet(:,k),cal.Zpos);
    zFocus(k) = out(2);
end
zFocus = zFocus(cal.neworder); %give the right order to the channels
distBetweenCamPlanes = abs(mean(diff(zFocus(1:2:end))) + mean(diff(zFocus(2:2:end))))/2;
target    = distBetweenCamPlanes/2;
distBetweenPlane = diff(zFocus);
offTarget1 = distBetweenPlane - target;
offTarget = mean(abs(offTarget1));

message = sprintf('The difference between the target and the current plane conformation is %d',offTarget);
disp(message);
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



