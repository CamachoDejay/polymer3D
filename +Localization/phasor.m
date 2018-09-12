function [row,col,e,MagX,MagY] = phasor(ROI)
%PHASOR calculates the sub-pixel localization of a single particle in the
%ROI. NOTE, x: position along the first dimention, y: position alog the
%second dimention.
%   Detailed explanation goes here

% I did a quick test and fft2 can work on integers and gives back doubles.
% assert(isa(ROI,'double'),'ROI image must be a double')

% size of ROI
roiSize = size(ROI);
% some tests about ROI
assert(roiSize(1)==roiSize(2), 'ROI must be a n-by-n matrix')
% This is not needed but for the moment I will leave it so it is compatible
% with gradient fit method.
assert(all(mod(roiSize,2)==[1,1]),'input ROI dimentions must be odd')

% Perform 2D fourier transformation on the ROI
fftROI = fft2(ROI);

% size of ROI
roiSize = roiSize(1);
% ROI center pixel
roiCenter = median([1,roiSize]);

% calculate the angle of the x-phasor from the first fourier coefficient in
% X;
angX = angle(fftROI(2,1));
% correct the angle
if angX>0
    angX = angX -2*pi;
end
% normalize the angle by 2pi and the ROI size
PositionX = (abs(angX) / (2*pi/roiSize) + 1);


% calculate the angle of the x-phasor from the first fourier coefficient in
% X;
angY = angle(fftROI(1,2));
% correct the angle
if angY>0
    angY = angY -2*pi;
end
% normalize the angle by 2pi and the ROI size
PositionY = (abs(angY) / (2*pi/roiSize) + 1);

% magnitude of the X and Y phasors
MagX = abs(fftROI(1,2));
MagY = abs(fftROI(2,1));

elip = MagX/MagY;

% this shift is so it is consistent with our other localization methods,
% which give the localize position in relation to the center of the ROI
row = PositionX-roiCenter;
col = PositionY-roiCenter;
e = elip;

end

