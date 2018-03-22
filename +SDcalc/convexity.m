function [ vals, names ] = convexity( contour )
%UNTITLED7 Summary of this function goes here
%   Input - must be a 2xn matrix containing the [x,y] coordinates of n
%   points, please note that convexity, being a function of perimeter, can
%   not be defined for contours with less than 3 points.

names = {'Convexity'};
vals = cell(length(names),1);

if isempty(contour)
    return
else


    assert(size(contour,1)==2,['Input must be a 2xn matrix containing'...
                               ' the [x,y] coordinates of n points'])
    assert(size(contour,2)>3,['Input must be a 2xn matrix containing'...
                               'the [x,y] coordinates of n points and a'...
                        ' perimeter cant be defined for less than 3 points'])

    k = convhull(contour(1,:),contour(2,:)); %k are indices
    CH = contour(:,k); % CH contains the coordinates of the convex hull

    [perim_shape]      = perimeter_li( contour );
    [perim_convexhull] = perimeter_li( CH );

    convex = perim_convexhull{1} / perim_shape{1};

    vals(1) = {convex};
end

end

