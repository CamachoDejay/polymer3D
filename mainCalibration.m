% this script is a work in progress, the plan is to fix here the quick
% callibration of the channels. Basically I want to be able to lead the
% movies and built a dataset of [xDim,yDim,channel,frame] where all
% channels represent the same region of the image to a pixel resolution.
clear
close all
clc

fPath = 'C:\Users\Boris\Dropbox\MultiPlane Data\2018\03 - Mar\ZstackCalibration-100nmBeads-640nm_2';
fName = 'ZstackCalibration-100nmBeads-640nm_2_MMStack_Pos0.ome.tif';

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
    [~,ind] = max(cal.focusMet(:,k));
    zFocus2(k) = cal.Zpos(ind);
end
dist2Firstplane = zFocus(1)-zFocus;
dist2Firstplane2 = zFocus2(1)-zFocus2;
distanceBetweenCamA = zFocus(1:4) - zFocus(5:8);
distanceBetweenCamB = zFocus(5:7) - zFocus(2:4);
distanceBetweenCam2A = zFocus2(1:4) - zFocus2(5:8);
distanceBetweenCam2B = zFocus2(5:7) - zFocus2(2:4);
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



