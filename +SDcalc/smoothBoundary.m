function [ boundary_S ] = smoothBoundary( boundary )
% Smoothing of boundaries as given by the function bwboundaries
%   Detailed explanation goes here

% check that boundaries come with expected format
assert(size(boundary,2) == 2, 'Boundary must consist of a 2 col matrix')
assert(size(boundary,1) >= 3, 'Boundary must consist of more than 2 points')
if boundary(1,:) == boundary(end,:)
    b = boundary(1:end-1,:);
%     disp('Working on boundary')
else
    warning('Boundary not given in expected format')
    boundary_S = [];
    return
end

% calculate number of boundary points
points_n_i = length(b);
% calculate how much to extend the array
margin     = round(points_n_i/10);

% extend boundary to avoid edge effects on smoothing
extend_1   =  b(end-margin+1:end,:);
extend_2   =  b(1:margin,:);
extended_b = [extend_1; b; extend_2];

% smooth the contour
extended_b(:,1) = smooth(extended_b(:,1));
extended_b(:,2) = smooth(extended_b(:,2));

% select only the correct part of the array (remove extension)
boundary_S = extended_b(margin+1:points_n_i+margin,:);
% close boundary to keep formating agreement
boundary_S = [boundary_S; boundary_S(1,:)];
end

