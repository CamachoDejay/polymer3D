function [ vals, names ] = solidity( contour )
% describes the extent to which the shape is convex or concave and it is
% deined by the ratio between the shape's area and that of its convexhull
%   Input - must be a 2xn matrix containing the [x,y] coordinates of n
%   points, please note that convexity, being a function of perimeter, can
%   not be defined for contours with less than 3 points.

names = {'Solidity'};
vals = cell(length(names),1);

if isempty(contour)
    return
else

    assert(size(contour,1)==2,['Input must be a 2xn matrix containing'...
                               ' the [x,y] coordinates of n points'])
    assert(size(contour,2)>3,['Input must be a 2xn matrix containing'...
                               'the [x,y] coordinates of n points and an'...
                        ' area cant be defined for less than 3 points'])

    k = convhull(contour(1,:),contour(2,:)); %k are indices
    CH = contour(:,k); % CH contains the coordinates of the convex hull

    [~,~,Area_shape]        = centroid_by_area( contour );
    [~,~,Area_convexhull]   = centroid_by_area( CH );

    Area_shape      = abs(Area_shape);
    Area_convexhull = abs(Area_convexhull);

    soli = Area_shape / Area_convexhull;

    vals(1) = {soli};

end
end

