% this script is a work in progress, the plan is to fix here the quick
% callibration of the channels. Basically I want to be able to lead the
% movies and built a dataset of [xDim,yDim,channel,frame] where all
% channels represent the same region of the image to a pixel resolution.
clear
close all
clc

% path to the callibration
fPath = '/Users/rafa/Documents/MATLAB/data/Boris/180322-Boris-Calmultiplane/BeadsCalibrationZStack_2';
fName = 'BeadsCalibrationZStack_2_MMStack_Pos0.ome.tif';

fPath = [fPath filesep fName];

% Calculate calibration
[cal] = mpSetup.cali.calculate(fPath, false);

disp('Done with calibration')
%%
% load and calibrate, when applied to the calibration data then we should
% be able to demonstrate that it works

fPath = '/Users/rafa/Documents/MATLAB/data/Boris/180322-Boris-Calmultiplane/TL-OD2-200msExposureR2_1';
fName = 'TL-OD2-200msExposureR2_1_MMStack_Pos0.ome.tif';
fPath = [fPath filesep fName];

ff = 1;
frame = 1:250;
[data, frameInfo, movInfo] = mpSetup.loadAndCal( fPath, cal, frame);

%%
% for ii = 1:250
%     for jj = 1:8
%         subplot(2,4,jj)
%         imagesc(data(:,:,jj,ii))
%         axis image
%     end
%     waitforbuttonpress
%     
% end

%% detect particles
% for a water immersion obj the best-fit gasuss to the PSF has 
% sigma = 0.25 wavelength / NA
objNA  = 1.2;
emWave = 600;
sigma_nm = 0.25 * emWave/objNA;
FWHMnm = sigma_nm * sqrt(8*log(2));         

GLRTprops.delta  = 6;
GLRTprops.pxSnm  = 100;
GLRTprops.FWHMnm = FWHMnm;
GLRTprops.chi2   = 80;

rTh = 5; % in pixels
ROIrad = 10;

pList = [];
p = [];

for fIdx = frame
    
    sfData = data(:,:,:,fIdx);
    imSize = size(sfData);
    % detect particles
    [consLoc,totLoc] = mpSetup.localize(sfData, rTh, GLRTprops);
    % build ROIs
    [ROIs] = Misc.getROIs(consLoc,ROIrad,imSize);
    
    for i = 1:size(consLoc,1)
        tmpLoc = consLoc(i,:);
        tmpROI = ROIs(i,:);
        p = mpSetup.particle(tmpLoc,fIdx,tmpROI,sfData);

        pList = [pList, p];
    end
    disp(['done for frame ' num2str(fIdx)])

end

%%
objNA    = 1.2;
emWave   = 600;
pxSizeNm = 100;
pList.setPSFprops(objNA, emWave, pxSizeNm);
%%

tic
pList.superResolve;
disp('Done with SR-loc')
toc
tic
p.superResolve;
toc

%%
test = findobj(pList,'frame',1);
tmpVal = cat(1,test.superResLoc);
tmpX = tmpVal(:,1);
tmpY = tmpVal(:,2);
tmpZ = tmpVal(:,3);
% scatter3 (tmpX,tmpY, tmpZ)
scatter (tmpX,tmpY,'kx')

% axis image
shg

cols = {'k','r','g','b','y'};
hold on
for ii = 2:250
    test = findobj(pList,'frame',ii);
    cidx = ceil(ii/50);
    tmpVal = cat(1,test.superResLoc);
    tmpX = tmpVal(:,1);
    tmpY = tmpVal(:,2);
    tmpZ = tmpVal(:,3);
%     scatter3 (tmpX,tmpY, tmpZ)
    scatter (tmpX,tmpY,[cols{cidx} 'x '])
    
end
hold off
axis image


% make a common list of sm detections? should I have a test for seeing a
% molecule in at least 3 (or X) planes? once I have the common list I have
% to cropt the ROIs [dx, dy, 8] and do the fine localization. We will for
% the moment pick the best as I do not have a super-resolved registration
% matrix between all channels. This should be generated in order to use
% info from multiple planes to increase fit accuracy. it is from z tacks of
% beads that we can create such a registration.



