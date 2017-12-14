function [ pos, inten ] = smDetection( im_in, delta, FWHM_pix, chi2 )
%smDetection finds the positions where a single molecule is most probably
%located in the input image using the generalized likelihood ratio test.
%   These postions should be further tested to ensure that a SM is in fact
%   there
%   See also GLRTFILTERING

% use GLRT to generate a new image where sharp peaks indicate the presence
% of a SM
BW = Localization.GLRTfiltering ( im_in, delta, FWHM_pix, chi2  );
% removing small objects - noise
BW = bwareaopen(BW, 4,8);
% find all isolated areas in the SM detection image
L  = bwlabel(BW);
% get the mean intensity and centroid of each detected molecule
stats = regionprops(L,im_in,'MeanIntensity','WeightedCentroid');

% generate outputs
pos      = cat(1,stats.WeightedCentroid);
inten    = cat(1,stats.MeanIntensity);

end

