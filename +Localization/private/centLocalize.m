function [x, y, e] = centLocalize(ROI)
% CENTLOCALIZE: estimates localization of the emitters psf [y,y] and its
% ellipticity [e] using centroid based method. ROI is a m-by-n matrix where
% m=n, and m,n must be odd;
%   Inputs:
%       ROI: Region of interest containing the centered image of a single
%       emitter.
%   Outputs:
%       x,y: Estimated position of the emitter.
%       e:   Estimated ellipticity of the emitters PSF.  
%   Author: Rafael Camacho camachodejay@yahoo.com, github user:
%   camachodejay. The centroid method idea comes from the publication
%   quickPALM and I borrowed some of the implementation from the work of
%   Hongqiang Ma (2015) which can be found in http://www.pitt.edu/~liuy.

% ROI dimensions: In imaging we ussually refer to the colums as X and to
% the rows as Y.
[sizeY, sizeX] = size(ROI);
% this function was designed for ROI with same x and y dimention and must
% be odd
assert(sizeY==sizeX,'ROI must be a m-by-n matrix where m=n')
assert(mod(sizeY,2)==1,'input ROI dimentions must be odd')

% radius of the ROI, we assume that the emitter's PSF is centered in the
% ROI. Thus we are interested in the 'radious of the ROI
% ROIrad = (sizeY-1)/2;

% pixIdx: pixel indices, notice that as x and y dim are the same I just
% have to define one vector for both X and Y
pixIdx = 1:sizeX;
% pixIdxY = 1:sizeY;

% Center pixel of ROI
centPix = median(pixIdx);

% centroid calculation
% ROI integral over the fisrt dimension, this means that the intensity
% trace obtained sROI_x is a function of x
sROI_x = sum(ROI,1);
% ROI integral over the second dimension, this means that the intensity
% trace obtained sROI_y is a function of y
sROI_y = sum(ROI,2);
% remove off set so that the min integral value is 0
% WARNING this step is detrimental to the ellipticity estimation if the PSF
% did not to a large extent fit into the ROI.
% sROI_x = sROI_x - min(sROI_x);
% sROI_y = sROI_y - min(sROI_y);
bgEstimate = min([sROI_x, sROI_y']);
sROI_x = sROI_x - bgEstimate;
sROI_y = sROI_y - bgEstimate;
% Centroid of an arbitrady 1D function
Xc = sum(sROI_x.*pixIdx)/sum(sROI_x);
Yc = sum(sROI_y.*pixIdx')/sum(sROI_y);

% Create a vector of the pixel indices shifted acording to the centroid,
% this way the centroid is now at (0,0)
shiftX = pixIdx-Xc;
shiftY = pixIdx-Yc;
% What is W? it seems to be a measurement of the width of the distribution.
% It could be realted to the second moment or quadrupole
Wx = sum(sROI_x.*abs(shiftX))/sum(sROI_x);
Wy = sum(sROI_y.*abs(shiftY'))/sum(sROI_y);

% they seem to be filtering some data out of the ROI integrals. They set
% the tails of the integrals to 0 according to the Ws
idx2zero = or(pixIdx<(Xc-3*Wx),pixIdx>(Xc+3*Wx));
sROI_x(idx2zero)=0;

idx2zero = or(pixIdx<(Yc-3*Wy),pixIdx>(Yc+3*Wy));
sROI_y(idx2zero)=0;

% refine centroid calculation after filtering
Xc = sum(sROI_x.*pixIdx)/sum(sROI_x);
Yc = sum(sROI_y.*pixIdx')/sum(sROI_y);
% Create a vector of the pixel indices shifted acording to the centroid,
% this way the centroid is now at (0,0)
shiftX = pixIdx-Xc;
shiftY = pixIdx-Yc;
% here we go again with the W
Wx = sum(sROI_x.*abs(shiftX))/sum(sROI_x);
Wy = sum(sROI_y.*abs(shiftY'))/sum(sROI_y);
% storing the centroid value asuming that the center of the ROI is 0,0
x = Xc-centPix;
y = Yc-centPix;
% e0 is the estimated ellipticity
e = Wy/Wx;

