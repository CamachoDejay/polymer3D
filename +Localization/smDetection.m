function [ pos, meanFAR, FAR ] = smDetection( im, delta, FWHM_pix, chi2 )
%smDetection finds the positions where a single molecule is most probably
%located in the input image using the generalized likelihood ratio test.
%   These postions should be further tested to ensure that a SM is in fact
%   there
%   See also GLRTFILTERING

% use GLRT to generate a new image where sharp peaks indicate the presence
% of a SM
% im2 = imgaussfilt(im,1);
[BW, ~, FAR] = GLRTfiltering ( im, delta, FWHM_pix, chi2  );
% removing small objects - noise
BW = bwareaopen(BW, 4,8);
% find all isolated areas in the SM detection image
L  = bwlabel(BW);
% get the mean intensity and centroid of each detected molecule
% stats = regionprops(L,im,'MeanIntensity','WeightedCentroid');
stats = regionprops(L,FAR,'MeanIntensity','WeightedCentroid');
% generate outputs
pos      = cat(1,stats.WeightedCentroid);
meanFAR    = cat(1,stats.MeanIntensity);

% note that due to the way that regionprops works the centroid positions
% are some what shifted eg if the image is of dimentions [200 500] you can
% end up with a pos at [500, 200]. This is due to MATLABs convention when
% working with images eg. imagesc. Moreover I want to keep things in such a
% way that the first pos index is the first im index. Thus I flip the col.

pos = flip(pos,2);
end

