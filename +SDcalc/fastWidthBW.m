function [ cell_width ] = fastWidthBW( im )
%CELLWIDTH Summary of this function goes here
%   c is the cell contour

% skelotanized image to look for center line
skeletonizedImage = bwmorph(im, 'skel', inf);
% eucledian distance from border image
edtImage          = bwdist(~im);
%now I find the values of the center line
distanceFromEdge = edtImage(skeletonizedImage);
% get those values that are larger than 0
d                = median(distanceFromEdge);
% calculate cell width
cell_width        = d*2; % stimate the width of the snake

end

