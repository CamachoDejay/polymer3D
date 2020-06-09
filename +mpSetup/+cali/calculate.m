function [cal, movInfo] = calculate(fPath,nChan, correctInt, flipCam2)
%CALCULATE calculates calibration for the multiplane setup. NOTE that if
%you choose to correct for intensity differences the data changes form
%uint16 to double because we have to multiply by a correction factor
%(double).

switch nargin
    case 1
        
        nChan = 4;
        flipCam2 = true;
        
    case 2
        % this is the normal way, I cant see how it would be different
        correctInt = true;
        flipCam2 = true;
        
    case 3
        
        flipCam2 = true;
        
end

cal.correctInt = correctInt;


cal.reorder    = true;


% load general information about the calibration movie
[frameInfo, movInfo, ~ ]= Load.Movie.ome.getInfo( fPath );

% load the data 
[ movC1, movC2, idx ] = Load.Movie.ome.load( frameInfo, movInfo, 1:movInfo.maxFrame );

h = waitbar(0,'Please wait...');
% store information about position given by the motors
posInfo = cat(1,frameInfo.Pos);
Z1 = posInfo(idx(:,1),3);
Z2 = posInfo(idx(:,2),3);

% if need be we flip camera 2, this is generally the case
cal.flipCam2 = flipCam2;
if cal.flipCam2
    movC2 = flip(movC2,2);
end

% % max projection image
% max_im1 = max(movC1,[],3);
% max_im2 = max(movC2,[],3);

meanIm1 = mean(movC1,3);
meanIm2 = mean(movC2,3);

waitbar(.1,h,'Finding channels')
% find channels
[ chCentCam1, ~, commonW1 ] = mpSetup.cali.findChannels( meanIm1, false, nChan );
[ chCentCam2, ~, commonW2 ] = mpSetup.cali.findChannels( meanIm2, false,nChan );

waitbar(.2,h,'getting ROIs')
% get ROI
commonwin = min([commonW1; commonW2]);
imS = size(meanIm1);
[ cal.ROI ] = mpSetup.cali.defineROI( commonwin, chCentCam1, chCentCam2, imS );

waitbar(.3,h,'getting channel data')
% get data for each channel identified. chData has dim im_size1 im_size2
% 4channels Nframes
[ chData1c, chData2c ] = mpSetup.cali.getChData( movC1, movC2, cal.ROI );

waitbar(.4,h,'getting focus metric')
% getting the focus metric as used in EPFL we might want to change this to
% a gradient method.
[ cal.focusMet, cal.inFocus, cal.fit ] = mpSetup.cali.getFocusMetric( chData1c, chData2c , Z1, Z2 );
cal.Zpos = Z1;

% figure(2)
% for i = 1:4
%     plot(Z1,focusMet(:,i))
%     hold on
% end
% for i = 5:8
%     plot(Z2,focusMet(:,i))
% end
% hold off

waitbar(.5,h,'getting new order for channels')
% find the new order for the camera channels
[ cal.neworder, cal.inFocus ] = mpSetup.cali.getNewOrder( cal.inFocus );

waitbar(.7,h,'getting image shifts')
% find image shift in order to have the same ROIs to a pixel resoltuon
[ imShifts ] = mpSetup.cali.simpleImShift2( cal.inFocus, chData1c, chData2c );

waitbar(.8,h,'refining ROIs')
% refine the ROIs to consider the shifts
[ cal.ROI ] = mpSetup.cali.refineROI( cal.ROI, imShifts );

figure()
subplot(2,1,1)
imagesc(meanIm1)
axis image
for i = 1:nChan
    rectangle('Position', cal.ROI(i,:))
end
title('Camera 1 with ROIs')
subplot(2,1,2)
imagesc(meanIm2)
axis image
for i = nChan+1:2*nChan
    rectangle('Position', cal.ROI(i,:))
end
title('Camera 2 with ROIs')

mpSetup.cali.plotCal(meanIm1,meanIm2, cal.ROI);

if cal.correctInt
    waitbar(.9,h,'Correcting intensity')
    % update the channel data
    [ chData1c, chData2c ] = mpSetup.cali.getChData( movC1, movC2, cal.ROI );
    % calculate intensity correction
    [ cal.Icorrf ] = mpSetup.cali.findChInt( chData1c, chData2c, cal.inFocus );
    maxInt = max(cal.fit(:,2:2:end),[],1);
    cal.Icorrf = maxInt./max(maxInt);
else
    cal.Icorrf = ones(8,1);
end

close(h)
end


