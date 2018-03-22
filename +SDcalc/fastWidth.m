function [ cell_width ] = fastWidth( c )
%CELLWIDTH Summary of this function goes here
%   c is the cell contour

assert(size(c,1) == 2, 'contour must be a [2xn] matrix' );
assert(size(c,2) > 2, 'a contour must constist of more than 2 points' );

x = c(1,:);
y = c(2,:);

x = x - min(x)+1;
y = y - min(y)+1;

x = round(x);
y = round(y);



im = false(max(y),max(x));

linearInd = sub2ind(size(im), y, x);
im(linearInd) = true;
im = imfill(im,'holes');
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

